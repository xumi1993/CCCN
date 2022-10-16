import obspy
import numpy as np
import os
from os.path import join
import pyasdf

def stack_all(path):
    ds = pyasdf.ASDFDataSet(path)
    sta = ds.waveforms.list()[0]
    tags = [tag for tag in ds.waveforms[sta].get_waveform_tags()]
    st_all = obspy.Stream()
    for tag in tags:
        st_all += ds.waveforms[sta][tag]
    st_all.stack()
    st_all.normalize()
    st_all[0].stats.starttime =  ds.waveforms[sta][tags[0]][0].stats.starttime
    return st_all

def symmetrize(st_all):
    for i in range(len(st_all)):
        mid_pos = int((st_all[i].stats.npts-1)/2)
        tr_sym = np.zeros(mid_pos+1)
        tr_sym[0] = st_all[i].data[mid_pos]
        tr_sym[1::] = st_all[i].data[0:mid_pos][::-1]+st_all[i].data[(mid_pos+1)::]
        st_all[i].data = tr_sym
