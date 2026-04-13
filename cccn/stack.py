import numpy as np
import obspy
from scipy.signal import hilbert
from scipy.fftpack import next_fast_len


def linear_stack(streams:obspy.core.stream.Stream):
    """
    Stack arrays of seismic traces by averaging them together. This method assumes that the traces are aligned in time and have the same sampling rate.
    """

    st = streams.copy()
    st.stack('linear')
    return st[0]


def phase_weighted_stack(streams:obspy.core.stream.Stream, power:int=2):
    """
    Phase-weighted stack (PWS) is a method of stacking seismic traces that emphasizes the coherent energy across multiple traces while suppressing incoherent noise. The PWS is calculated by first computing the analytic signal of each trace using the Hilbert transform, then extracting the instantaneous phase and amplitude. The phase information is used to weight the amplitude of each trace, giving more weight to traces that are in phase with each other.

    Parameters:
    streams (obspy.core.stream.Stream): A stream of seismic traces to be stacked.
    power (int): The exponent used to weight the phase information. Higher values will give more weight to traces that are in phase.

    Returns:
    np.ndarray: The phase-weighted stack of the input traces.
    """
    arr = np.array([tr.data for tr in streams])
    assert arr.ndim == 2, "Traces must have equal length for stacking"

    # Compute the analytic signal for each trace
    analytic_signals = hilbert(arr, N=next_fast_len(arr.shape[1]), axis=1)
    
    # Extract the instantaneous phase and amplitude
    phase = np.angle(analytic_signals)
    phase_stack = np.abs(np.mean(np.exp(1j * phase), axis=0)) ** power
    
    # Compute the weighted stack
    weighted = np.multiply(arr, phase_stack)

    stacked = streams[0].copy()
    stacked.data = np.mean(weighted, axis=0)

    return stacked


def selective_stack(streams:obspy.core.stream.Stream, epsilon:float=1e-3):
    """
    Selective stacking is a method of stacking seismic traces that involves selecting a subset of traces based on certain criteria, such as signal-to-noise ratio or coherence. The selected traces are then stacked together to enhance the signal while suppressing noise.

    Parameters:
    streams (obspy.core.stream.Stream): A stream of seismic traces to be stacked.
    epsilon (float): The threshold for selecting traces, defaults to 1e-3.

    Returns:
    trace: The selectively stacked trace.
    """
    arr = np.array([tr.data for tr in streams])
    assert arr.ndim == 2, "Traces must have equal length for stacking"

    cc = np.ones(arr.shape[0])
    newstack = np.mean(arr, axis=0)
    for i, trace in enumerate(arr):
        cc[i] = np.sum(np.multiply(newstack, trace.T))
    ik = np.where(cc >= epsilon)[0]
    newstack = np.mean(arr[ik, :], axis=0)
    stacked = streams[0].copy()
    stacked.data = newstack
    return stacked