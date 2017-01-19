import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import obspy
from obspy.signal.util import next_pow_2
from make_one_folder import smooth


def whiten(data, nt, Nfft, delta, freqmin, freqmax, plot=False):
    """This function takes 1-dimensional *data* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *freqmin* and *freqmax*
    and returns the whitened fft.

    :type data: :class:`numpy.ndarray`
    :param data: Contains the 1D time series to whiten
    :type Nfft: int
    :param Nfft: The number of points to compute the FFT
    :type delta: float
    :param delta: The sampling frequency of the `data`
    :type freqmin: float
    :param freqmin: The lower frequency bound
    :type freqmax: float
    :param freqmax: The upper frequency bound
    :type plot: bool
    :param plot: Whether to show a raw plot of the action (default: False)

    :rtype: :class:`numpy.ndarray`
    :returns: The FFT of the input trace, whitened between the frequency bounds
"""

    if plot:
        fig = plt.figure(figsize=(20, 32), dpi=60)
        plt.subplot(411)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    Napod = 100
    Nfft = int(Nfft)
    dom = 1/delta/Nfft
    freqVec = scipy.fftpack.fftfreq(Nfft,d=delta)[:Nfft/2]
    
    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    low = J[0] - Napod
    if low <= 0:
        low = 1

    porte1 = J[0]
    porte2 = J[-1]
    high = J[-1] + Napod
    if high > Nfft / 2:
        high = Nfft // 2

    FFTRawSign = scipy.fftpack.fft(data, Nfft)
    
    if plot:
        plt.subplot(412)
        axis = np.arange(len(FFTRawSign)) * dom
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')

    # Left tapering:
    FFTRawSign[0:low] *= 0
    FFTRawSign[low:porte1] = np.cos(np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(1j * np.angle(FFTRawSign[low:porte1]))
    # Pass band:
    FFTRawSign[porte1:porte2] = np.exp(1j * np.angle(FFTRawSign[porte1:porte2]))
    # Right tapering:
    FFTRawSign[porte2:high] = np.cos(np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(1j * np.angle(FFTRawSign[porte2:high]))
    FFTRawSign[high:Nfft+1] *= 0
    
    # Hermitian symmetry (because the input is real)
    FFTRawSign[-Nfft/2+1:] = FFTRawSign[1:Nfft/2].conjugate()[::-1]

    if plot:
        plt.subplot(413)
        axis = np.arange(len(FFTRawSign)) * dom
        plt.axvline(low, c='g')
        plt.axvline(porte1, c='g')
        plt.axvline(porte2, c='r')
        plt.axvline(high, c='r')

        plt.axvline(Nfft - high, c='r')
        plt.axvline(Nfft - porte2, c='r')
        plt.axvline(Nfft - porte1, c='g')
        plt.axvline(Nfft - low, c='g')

        plt.plot(axis, np.abs(FFTRawSign))
        plt.xlim(0, max(axis))
        plt.ylim(0, 1.5)
        plt.title("Amp after Whiten")

        wdata = np.real(scipy.fftpack.ifft(FFTRawSign))
        plt.subplot(414)
        plt.plot(np.arange(len(wdata)) * delta, wdata)
        plt.xlim(0, nt * delta)
        plt.title("Output trace after Whiten and filter")
        plt.show()
    
    return FFTRawSign


if __name__ == '__main__':
    import time
    a = obspy.read('IC.BJT.00.BHZ.norm')[0]
    nt = a.stats.npts
    N = next_pow_2(nt)
#    a = np.sin(a) + np.sin(a / 4.) + np.sin(a / 16.)
#    a -= a.mean()
    t = time.clock()
    for i in range(1):
        whiten(a.copy(), nt, N, a.stats.delta, 0.02, 0.067, plot=False)
    print("1000 loops:", (time.clock()-t) * 1000, "ms")
    fftdata = whiten(a.copy(), nt, N, a.stats.delta, 0.02, 0.067, plot=True)
'''
    wdata = np.real(scipy.fftpack.ifft(fftdata))
    wdata = wdata[0:nt]/smooth(np.abs(wdata[0:nt]), half_len=20)
    plt.plot(np.arange(len(wdata)) * a.stats.delta, wdata)
    plt.xlim(0, nt*a.stats.delta)
    plt.show()
'''



