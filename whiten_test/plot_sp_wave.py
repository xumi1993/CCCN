#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import obspy

fig = plt.figure(figsize=(20, 8), dpi=60)

sp = obspy.read("out.[ri][lm]")
sp_amp = np.abs(sp[0].data+1j*sp[1].data)
axis = np.arange(sp[0].stats.npts) * sp[0].stats.delta
plt.subplot(211)
plt.plot(axis, sp_amp)
plt.xlim(0, max(axis))
plt.title("spectrum after whiten (filter4)")

tr = obspy.read("out.tr")[0]
axis = np.arange(tr.stats.npts) * tr.stats.delta
plt.subplot(212)
plt.plot(axis, tr.data)
plt.xlim(0, max(axis))
plt.title("time trace after whiten (filter4)")

plt.show()
