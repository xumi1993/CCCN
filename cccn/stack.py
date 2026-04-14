import numpy as np
import obspy
from obspy.io.sac import SACTrace
from scipy.signal import hilbert
from scipy.fftpack import next_fast_len
import h5py


def _as_2d_array(arrays):
    arr = np.asarray(arrays)
    if arr.ndim != 2:
        raise ValueError('Input must be a 2D numpy array: (n_traces, n_samples)')
    return arr


def linear_stack(arrays):
    """
    Stack 2D arrays of seismic traces by averaging them together.
    """
    arr = _as_2d_array(arrays)
    return np.mean(arr, axis=0)


def phase_weighted_stack(arrays, power:int=2):
    """
    Phase-weighted stack (PWS) is a method of stacking seismic traces that emphasizes the coherent energy across multiple traces while suppressing incoherent noise. The PWS is calculated by first computing the analytic signal of each trace using the Hilbert transform, then extracting the instantaneous phase and amplitude. The phase information is used to weight the amplitude of each trace, giving more weight to traces that are in phase with each other.

    Parameters:
    arrays (np.ndarray): A 2D array of shape (n_traces, n_samples).
    power (int): The exponent used to weight the phase information. Higher values will give more weight to traces that are in phase.

    Returns:
    np.ndarray: The phase-weighted stack of the input traces.
    """
    arr = _as_2d_array(arrays)
    n_samples = arr.shape[1]

    # Compute the analytic signal for each trace
    analytic_signals = hilbert(arr, N=next_fast_len(n_samples), axis=1)
    
    # Extract the instantaneous phase and amplitude
    phase = np.angle(analytic_signals)
    phase_stack = np.abs(np.mean(np.exp(1j * phase), axis=0)) ** power
    # Hilbert with N>n_samples pads output; trim back to original trace length.
    phase_stack = phase_stack[:n_samples]
    
    # Compute the weighted stack
    weighted = np.multiply(arr, phase_stack)
    return np.mean(weighted, axis=0)


def selective_stack(arrays, epsilon:float=1e-3):
    """
    Selective stacking is a method of stacking seismic traces that involves selecting a subset of traces based on certain criteria, such as signal-to-noise ratio or coherence. The selected traces are then stacked together to enhance the signal while suppressing noise.

    Parameters:
    arrays (np.ndarray): A 2D array of shape (n_traces, n_samples).
    epsilon (float): The threshold for selecting traces, defaults to 1e-3.

    Returns:
    np.ndarray: The selectively stacked trace.
    """
    arr = _as_2d_array(arrays)

    cc = np.ones(arr.shape[0])
    newstack = np.mean(arr, axis=0)
    for i, trace in enumerate(arr):
        cc[i] = np.sum(np.multiply(newstack, trace.T))
    ik = np.where(cc >= epsilon)[0]
    if len(ik) == 0:
        return newstack
    return np.mean(arr[ik, :], axis=0)


def read_cor(fname: str, sort_tags: bool = True):
    """Read a COR_*.h5 file and return header, tags and stacked-ready arrays."""
    with h5py.File(fname, 'r') as f:
        header = dict(f.attrs)
        tags = list(f.keys())
        if sort_tags:
            tags = sorted(tags)
        arrays = [f[tag][:] for tag in tags]

    if arrays:
        data = np.asarray(arrays)
    else:
        data = np.empty((0, 0))

    return {
        'header': header,
        'tags': tags,
        'arrays': data,
    }


class Stack:
    def __init__(self, fname: str, sort_tags: bool = True) -> None:
        self.fname = fname
        data = read_cor(self.fname, sort_tags=sort_tags)
        self.header = data['header']
        self.tags = data['tags']
        self.arrays = data['arrays']
        self.stacked = obspy.Trace()

    def _build_trace(self, data):
        data = np.asarray(data, dtype=np.float32)
        sactr = SACTrace(data=data, npts=len(data))
        sactr.delta = float(self.header.get('delta'))
        sactr.kevnm = f"{self.header.get('network1', '')}.{self.header.get('station1', '')}"
        sactr.kcmpnm = f"{self.header.get('channel1', '')}_{self.header.get('channel2', '')}"
        sactr.knetwk = self.header.get('network2', '')
        sactr.kstnm = self.header.get('station2', '')
        tr = sactr.to_obspy_trace()
        return tr

    def stack(self, method: str = 'linear', **kwargs):
        if self.arrays.size == 0:
            raise ValueError('No CCF data found in file')

        m = method.lower()
        if m == 'linear':
            stacked = linear_stack(self.arrays)
            self.stacked = self._build_trace(stacked)
        elif m == 'pws':
            power = kwargs.get('power', 2)
            stacked = phase_weighted_stack(self.arrays, power=power)
            self.stacked = self._build_trace(stacked)
        elif m == 'selective':
            epsilon = kwargs.get('epsilon', 1e-3)
            stacked = selective_stack(self.arrays, epsilon=epsilon)
            self.stacked = self._build_trace(stacked)
        else:
            raise ValueError("method must be one of: 'linear', 'pws', 'selective', got: {}".format(method))
    
    def fold(self):
        if (not self.stacked) or (self.stacked.data is None):
            raise ValueError('No stacked trace found. Please run stack() method first.')
        mid_pos = int((self.stacked.stats.npts-1)/2)
        sym = np.zeros(mid_pos+1)
        sym[0] = self.stacked.data[mid_pos]
        sym[1::] = self.stacked.data[0:mid_pos][::-1]+self.stacked.data[(mid_pos+1)::]
        tr_sym = self._build_trace(sym)
        return tr_sym