import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
import obspy
from obspy.signal.util import next_pow_2
from scipy.signal import correlate
from whiten import whiten


sacfile = "../example/2016.01.01/IC.*.00.BHZ.norm"
f1 = 0.01
f2 = 0.02
f3 = 0.067
f4 = 0.08
st = obspy.read(sacfile)
nt = st[0].stats.npts
nft = int(next_pow_2(nt))
dt = st[0].stats.delta
ftr1 = whiten(st[0].data,nft,dt,f1,f2,f3,f4)
ftr2 = whiten(st[1].data,nft,dt,f1,f2,f3,f4)
tr1 = fftpack.ifft(ftr1, nt)
tr2 = fftpack.ifft(ftr2, nt)

ccf = fftpack.ifft(ftr1*np.conj(ftr2), nt).real/nt
tcorr = np.arange(-nt + 1, nt)
dn = np.where(np.abs(tcorr) <= (nt-1)/2)[0]
#ccf = fftpack.ifftshift(ccf)
ccf = np.concatenate((ccf[-nt + 1:], ccf[:nt + 1]))
print(ccf.shape)
ccf = ccf[dn]
cci = correlate(tr1,tr2,"same")/nt
plt.subplot(211)
plt.plot(ccf)
plt.title("cc in freq domain")

plt.subplot(212)
plt.plot(cci)
plt.title("cc using scipy.signal.correlate")

plt.show()
