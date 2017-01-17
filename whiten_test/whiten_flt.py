import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import obspy
from obspy.signal.util import next_pow_2
from make_one_folder import smooth


def whiten(data, nt, Nfft, delta, f1, f2, f3, f4, plot=False):
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
    nt1 = f1/dom
    nt2 = f2/dom
    nt3 = f3/dom
    nt4 = f4/dom

    FFTRawSign = scipy.fftpack.fft(data, Nfft)
    
    if plot:
        plt.subplot(412)
        axis = np.arange(len(FFTRawSign)) * dom
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')
    
    FFTRawSign /= smooth(abs(FFTRawSign), half_len=20)
    # Left tapering:
    FFTRawSign[0:nt1] *= 0
    FFTRawSign[nt1:nt2] = np.cos(np.linspace(np.pi / 2., np.pi, nt2 - nt1+1)) ** 2 * FFTRawSign[nt1:nt2]
    # Pass band:
#    FFTRawSign[porte1:porte2] = np.exp(1j * np.angle(FFTRawSign[porte1:porte2]))
    # Right tapering:
    FFTRawSign[nt3:nt4] = np.cos(np.linspace(0., np.pi / 2., nt4 - nt3+1)) ** 2 * FFTRawSign[nt3:nt4]
    FFTRawSign[nt4:Nfft+1] *= 0
    
    # Hermitian symmetry (because the input is real)
    FFTRawSign[-Nfft/2+1:] = FFTRawSign[1:Nfft/2].conjugate()[::-1]

    if plot:
        plt.subplot(413)
        axis = np.arange(len(FFTRawSign)) * dom
        plt.axvline(nt1, c='g')
        plt.axvline(nt2, c='g')
        plt.axvline(nt3, c='r')
        plt.axvline(nt4, c='r')

        plt.axvline(Nfft - nt1, c='r')
        plt.axvline(Nfft - nt2, c='r')
        plt.axvline(Nfft - nt3, c='g')
        plt.axvline(Nfft - nt4, c='g')

        plt.plot(axis, np.abs(FFTRawSign))
        plt.xlim(0, max(axis))
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
    fftdata = whiten(a.copy(), nt, N, a.stats.delta, 0.01, 0.02, 0.067, 0.08, plot=True)
    print("whiten the trace:", (time.clock()-t) * 1000, "ms")
'''
    wdata = np.real(scipy.fftpack.ifft(fftdata))
    wdata = wdata[0:nt]/smooth(np.abs(wdata[0:nt]), half_len=20)
    plt.plot(np.arange(len(wdata)) * a.stats.delta, wdata)
    plt.xlim(0, nt*a.stats.delta)
    plt.show()
'''



