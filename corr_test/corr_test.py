import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
from obspy.signal.util import next_pow_2
from scipy.signal import correlate

x = np.arange(0,100,0.1)
lag = 200
Nt = len(x)
ccc = np.zeros(2*lag+1)
a = np.sin(x)
b = np.cos(x)
nft = int(next_pow_2(len(x)))
af = fftpack.fft(a,nft)
bf = fftpack.fft(b,nft)
cc = fftpack.ifft(af*np.conj(bf),len(x)).real/Nt
cci = fftpack.ifftshift(cc)
cc = fftpack.ifftshift(cc)
tcorr = np.arange(-Nt + 1, Nt)
dN = np.where(np.abs(tcorr) <= lag)[0]
cc = cc[dN]


plt.subplot(311)
plt.plot(cci)

plt.subplot(312)
plt.plot(cc)

plt.subplot(313)
plt.plot(correlate(a,b))

plt.show()
