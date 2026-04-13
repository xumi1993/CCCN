import numpy as np
import obspy
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

    # Compute the analytic signal for each trace
    analytic_signals = hilbert(arr, N=next_fast_len(arr.shape[1]), axis=1)
    
    # Extract the instantaneous phase and amplitude
    phase = np.angle(analytic_signals)
    phase_stack = np.abs(np.mean(np.exp(1j * phase), axis=0)) ** power
    
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
        self.stream = self.to_stream()


    def _build_trace(self, data):
        tr = obspy.Trace(data=np.asarray(data, dtype=np.float64))
        tr.stats.delta = float(self.header.get('delta', 1.0))

        # Map station/channel metadata from COR header into Trace stats.
        tr.stats.network = str(self.header.get('network1', ''))
        tr.stats.station = str(self.header.get('station1', ''))
        tr.stats.channel = '{}{}'.format(
            self.header.get('channel1', ''),
            self.header.get('channel2', ''),
        )

        tr.stats.cc = {
            'network1': self.header.get('network1', ''),
            'station1': self.header.get('station1', ''),
            'channel1': self.header.get('channel1', ''),
            'network2': self.header.get('network2', ''),
            'station2': self.header.get('station2', ''),
            'channel2': self.header.get('channel2', ''),
            'lag': self.header.get('lag', None),
        }
        return tr

    def stack(self, method: str = 'linear', **kwargs):
        if self.arrays.size == 0:
            raise ValueError('No CCF data found in file')

        m = method.lower()
        if m == 'linear':
            stacked = linear_stack(self.arrays)
            return self._build_trace(stacked)
        if m == 'pws':
            power = kwargs.get('power', 2)
            stacked = phase_weighted_stack(self.arrays, power=power)
            return self._build_trace(stacked)
        if m == 'selective':
            epsilon = kwargs.get('epsilon', 1e-3)
            stacked = selective_stack(self.arrays, epsilon=epsilon)
            return self._build_trace(stacked)
        raise ValueError("method must be one of: 'linear', 'pws', 'selective'")
