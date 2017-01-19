import obspy
import numpy as np
import os
from os.path import join

def stack_all(path,sta_pair):
    st_all = obspy.Stream()
    for each_pair in sta_pair:
        st = obspy.read(join(path,"*/COR",each_pair))
        ns = len(st)
        npts = st[0].stats.npts
        data = np.zeros([ns, npts])
        for i in range(ns):
            data[i,:] = st[i].data
        tmp_tr = st[0].copy()
        tmp_tr.data = data.sum(axis=0)
        st_all.append(tmp_tr)
    return st_all
