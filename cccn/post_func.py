import obspy
import numpy as np
import os
from os.path import join

def stack_all(path,sta_pair):
    st_all = obspy.Stream()
    for each_pair in sta_pair:
        st = obspy.read(join(path,"*/COR","COR_"+each_pair+".SAC"))
        ns = len(st)
        npts = st[0].stats.npts
        data = np.zeros([ns, npts])
        for i in range(ns):
            data[i,:] = st[i].data
        tmp_tr = st[0].copy()
        tmp_tr.data = data.sum(axis=0)
        st_all.append(tmp_tr)
    return st_all

def symmetrize(st_all):
    for i in range(len(st_all)):
        mid_pos = int((st_all[i].stats.npts-1)/2)
        tr_sym = np.zeros(mid_pos+1)
        tr_sym[0] = st_all[i].data[mid_pos]
        tr_sym[1::] = st_all[i].data[0:mid_pos][::-1]+st_all[i].data[(mid_pos+1)::]
        st_all[i].data = tr_sym
