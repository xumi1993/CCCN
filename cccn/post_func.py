import obspy
import numpy as np
from obspy.io.sac import SACTrace
from obspy.core.util.attribdict import AttribDict
from os.path import join
import pyasdf


class PostProcForNoise:
    def __init__(self, fname) -> None:
        self.fname = fname
        self.ds = pyasdf.ASDFDataSet(self.fname)
        self.stname = self.ds.waveforms.list()[0]
        self.tags = [tag for tag in self.ds.waveforms[self.stname].get_waveform_tags()]
        self.ccfstream = obspy.Stream()
        for tag in self.tags:
            self.ccfstream += self.ds.waveforms[self.stname][tag]
        self.evname = self.ds.events[0].resource_id.id.split('/')[-1]
        stainfo = self.ds.get_all_coordinates()
        self.evla = self.ds.events[0].origins[0].latitude
        self.evlo = self.ds.events[0].origins[0].longitude
        self.stla = stainfo[self.stname]['latitude']
        self.stlo = stainfo[self.stname]['longitude']

    def stack_all(self, normalize=True, add_sac_header=True):
        stack_st = self.ccfstream.copy().stack()
        if normalize:
            stack_st.normalize()
        # st_all[0].stats.update(ds.waveforms[sta][tag][0].stats.sac)
        stack_st[0].stats.starttime = self.ds.waveforms[self.stname][self.tags[0]][0].stats.starttime
        if add_sac_header:
            sacheader = SACTrace()._header
            sacheader['knetwk'] = self.stname.split('.')[0]
            sacheader['kstnm'] = self.stname.split('.')[1]
            sacheader['stla'] = self.stla
            sacheader['stlo'] = self.stlo
            sacheader['kevnm'] = self.evname
            sacheader['evla'] = self.evla
            sacheader['evlo'] = self.evlo
            stack_st[0].stats.sac = AttribDict(sacheader)
        return stack_st

    def symmetrize(self):
        for i in range(len(self.ccfstream)):
            mid_pos = int((self.ccfstream[i].stats.npts-1)/2)
            tr_sym = np.zeros(mid_pos+1)
            tr_sym[0] = self.ccfstream[i].data[mid_pos]
            tr_sym[1::] = self.ccfstream[i].data[0:mid_pos][::-1]+self.ccfstream[i].data[(mid_pos+1)::]
            self.ccfstream[i].data = tr_sym
