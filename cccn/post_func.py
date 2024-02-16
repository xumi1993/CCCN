import obspy
import numpy as np
from obspy.io.sac import SACTrace
from obspy.core.util.attribdict import AttribDict
from obspy.geodetics import gps2dist_azimuth
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
        self.evel = -self.ds.events[0].origins[0].depth
        self.stla = stainfo[self.stname]['latitude']
        self.stlo = stainfo[self.stname]['longitude']
        self.stel = stainfo[self.stname]['elevation_in_m']

    def stack_all(self, normalize=True, add_sac_header=True):
        self.stack_st = self.ccfstream.copy().stack()
        if normalize:
            self.stack_st.normalize()
        # st_all[0].stats.update(ds.waveforms[sta][tag][0].stats.sac)
        self.stack_st[0].stats.starttime = self.ds.waveforms[self.stname][self.tags[0]][0].stats.starttime
        if add_sac_header:
            sacheader = SACTrace()._header
            sacheader['knetwk'] = self.stname.split('.')[0]
            sacheader['kstnm'] = self.stname.split('.')[1]
            sacheader['stla'] = self.stla
            sacheader['stlo'] = self.stlo
            sacheader['kevnm'] = self.evname
            sacheader['evla'] = self.evla
            sacheader['evlo'] = self.evlo
            dis, _, _ = gps2dist_azimuth(self.evla, self.evlo, self.stla, self.stlo)
            sacheader['dist'] = dis/1000
            self.stack_st[0].stats.sac = AttribDict(sacheader)
        return self.stack_st

    def write_stack(self, outpath='./', format='txt'):
        """Write stacked CCF to SAC file or ASCII file

        Parameters
        ----------
        outpath : str, optional
            Path to saving, by default './'
        format : str, optional
            Format of files, ``sac`` or ``txt`` are available, by default ``txt``
        """
        if not hasattr(self, 'stack_st'):
            raise ValueError('Please stack CCF first')
        fname = '{}_{}-{}_{}'.format(
                self.stack_st[0].stats.channel,
                self.evname, self.stname, len(self.ccfstream))
        if format.lower() == 'sac':
            self.stack_st[0].write(join(outpath, '{}.sac'.format(fname)))
        elif format.lower() == 'txt':
            mid_pos = int((self.stack_st[0].stats.npts-1)/2)
            times = np.arange(mid_pos)*self.stack_st[0].stats.delta
            ccf1 = np.flip(self.stack_st[0].data[0:mid_pos])
            ccf2 = self.stack_st[0].data[mid_pos+1:]
            zeroamp = (self.stack_st[0].data[mid_pos-1:mid_pos+2])/3
            ccf1[0] = zeroamp
            ccf2[0] = zeroamp
            data = np.zeros([mid_pos+2, 3])
            data[0] = self.stlo, self.stla, self.stel
            data[1] = self.evlo, self.evla, self.evel
            data[2:, 0] = times
            data[2:, 1] = ccf1
            data[2:, 2] = ccf2
            np.savetxt(join(outpath, '{}.dat'.format(fname)), data)
        else:
            pass

    def symmetrize(self):
        for i in range(len(self.ccfstream)):
            mid_pos = int((self.ccfstream[i].stats.npts-1)/2)
            tr_sym = np.zeros(mid_pos+1)
            tr_sym[0] = self.ccfstream[i].data[mid_pos]
            tr_sym[1::] = self.ccfstream[i].data[0:mid_pos][::-1]+self.ccfstream[i].data[(mid_pos+1)::]
            self.ccfstream[i].data = tr_sym
