import obspy
import os
from os.path import join, exists
import numpy as np
from obspy.signal.util import next_pow_2
from obspy.core.inventory import Inventory, Network, Station, Channel
from scipy import fftpack
from obspy.core.event.event import Event
from obspy.core.event.origin import Origin
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.io.sac import SACTrace
from .util import *
from .para import Para
import time
import pyasdf
import glob
from .parallel import MyMPI
from .logger import Logger
from .io import NoiseFile


class CrossCorrelation():
    def __init__(self, para=None, level='INFO') -> None:
        self.mpi = MyMPI()
        self.log = Logger('CCCN', level=level)
        self.logger = self.log.get_logger()
        if para is None:
            self.para = Para()
        else:
            self.para = para
        self.rawst = obspy.Stream()
        self.para.bcast(self.mpi)
        self.dt = None
        self.nft = None
        self.ntr = None
        self.reftime = None
    
    def read_data(self, matchstr='*'):
        if self.mpi.world_rank == 0:
            self.logger.info(f'Reading {self.para.suffix} files in {self.para.datapath}')
            if exists(self.para.datapath):
                fname = join(self.para.datapath,f'{matchstr}*.{self.para.suffix}')
            else:
                fname = self.para.datapath
            try:
                self.rawst = obspy.read(fname)
            except:
                return False
            self.erease_duplicates()
            if self.para.target_dt is None:
                self.dt = self.rawst[0].stats.delta
            else:
                self.dt = self.para.target_dt
                if self.dt > self.rawst[0].stats.delta:
                    self.rawst.filter('lowpass', freq=1/(2*self.dt), corners=2, zerophase=True)
                self.rawst.resample(1/self.para.target_dt, no_filter=True)
        self.mpi.synchronize_all()
        self.dt = self.mpi.bcast(self.dt)
        return True

    def prepare_stainfo(self):
        self.network = []
        self.station = []
        self.channel = []
        self.staid = []
        if self.mpi.world_rank == 0:
            for i, tr in enumerate(self.rawst):
                self.network.append(tr.stats.network)
                self.station.append(tr.stats.station)
                self.channel.append(tr.stats.channel)
                self.staid.append(tr.stats.network+'_'+tr.stats.station+'_'+tr.stats.channel[-1])
        self.network = self.mpi.bcast(self.network)
        self.station = self.mpi.bcast(self.station)
        self.channel = self.mpi.bcast(self.channel)
        self.staid = self.mpi.bcast(self.staid)

    def erease_duplicates(self):
        if self.mpi.world_rank == 0:
            ids = [tr.stats.network+'_'+tr.stats.station+'_'+tr.stats.channel for tr in self.rawst]
            duplicated = list(set([x for x in ids if ids.count(x) > 1]))
            self.rawst = obspy.Stream([tr for tr in self.rawst if tr.stats.network+'_'+tr.stats.station+'_'+tr.stats.channel not in duplicated])

    def perwhiten(self):
        if self.mpi.world_rank == 0:
            if not self.rawst:
                return
            if self.para.reftime == 'day':
                self.reftime = obspy.UTCDateTime(self.rawst[0].stats.starttime.date)
            elif self.para.reftime == 'hour':
                self.reftime = obspy.UTCDateTime(self.rawst[0].stats.starttime.strftime('%Y%m%d%H'))
            elif self.para.reftime == 'minute':
                self.reftime = obspy.UTCDateTime(self.rawst[0].stats.starttime.strftime('%Y%m%d%H%M'))
            timestart = self.para.timeduration*self.para.cut_precentatge
            timeend =  self.para.timeduration*(1-self.para.cut_precentatge)
            self.nft = int(next_pow_2((timeend-timestart)/self.dt))
            fftst = self.rawst.copy()
            self.logger.info(f'Pre-processing of {self.reftime}')
            for tr in fftst:
                self.logger.debug(f'Pre-processing {tr.stats.network}_{tr.stats.station}_{tr.stats.channel}')
                #------- cut waveform ------ 
                cutbtime = self.reftime+timestart
                cutetime = self.reftime+timeend
                if tr.stats.starttime > cutbtime or tr.stats.endtime < cutetime:
                    fftst.remove(tr)
                    self.rawst.remove(tr)
                    continue
                tr.trim(cutbtime, cutetime)
                if np.isnan(tr.data).any() or tr.data.all() == 0:
                    self.logger.warning(f'Invalid value in {tr.stats.network}_{tr.stats.station}_{tr.stats.channel}, removed')
                    fftst.remove(tr)
                    self.rawst.remove(tr)
                    continue
                tr.detrend('linear')
                tr.detrend('constant')
                #----------normalize----------
                if self.para.nsmooth == 0:
                    tr.data /= np.abs(tr.data)
                elif self.para.nsmooth > 0:
                    tr.data /= smooth(np.abs(tr.data),half_len=self.para.nsmooth)
                else:
                    raise ValueError("Half window length must be greater than zero")
                #----------- Whiten -----------
                f1 = 1/(1.5*(1/self.para.freqmin))
                f4 = 1/(0.75*(1/self.para.freqmax))
                (tr.data,tr.stats.delta) = whiten(tr.data, self.nft, self.dt, f1, self.para.freqmin, self.para.freqmax, f4)
            self.ntr = len(fftst)
        self.mpi.synchronize_all()
        self.ntr = self.mpi.bcast(self.ntr)
        self.nft = self.mpi.bcast(self.nft)
        self.reftime = self.mpi.bcast(self.reftime)

        self.df = 1/self.dt/self.nft
        self.fftarr, self.fft_win = self.mpi.prepare_shm([self.ntr, self.nft], np.complex128)
        if self.mpi.world_rank == 0:
            for i, tr in enumerate(fftst):
                self.fftarr[i] = tr.data
        self.mpi.sync_from_main(self.fftarr)
        self.prepare_stainfo()

    def prepare_cc(self, action='w'):
        self.idxij = []
        if self.mpi.world_rank == 0:
            if self.para.src_mask:
                for src in self.para.src_mask:
                    try:
                        i = self.staid.index(src)
                    except:
                        continue
                    for j in np.arange(self.ntr):
                        if i != j:
                            self.idxij.append([i, j])
                # remove duplicates i, j and j, i
                self.idxij = [list(x) for x in set([tuple(sorted(x)) for x in self.idxij])]
            else:
                for i in np.arange(self.ntr-1):
                    for j in np.arange(i+1, self.ntr):
                        self.idxij.append([i, j])
        self.mpi.synchronize_all()
        self.idxij = self.mpi.bcast(self.idxij)
        for i, idxij in enumerate(self.idxij):
            if i % self.mpi.world_size == self.mpi.world_rank:
                channel = self.channel[idxij[0]][-1]+self.channel[idxij[1]][-1]
                sta_pair = '{}_{}'.format(f'{self.network[idxij[0]]}_{self.station[idxij[0]]}',
                                          f'{self.network[idxij[1]]}_{self.station[idxij[1]]}')
                self.create_dataset('COR_{}_{}.h5'.format(sta_pair, channel), idxij, action=action)
        self.mpi.synchronize_all()

    def run_cc(self, action='w'):
        self.prepare_cc(action)
        mid_pos = int(self.nft/2)
        nlag = int(self.para.maxlag/self.dt)
        for i, idxij in enumerate(self.idxij):
            if i % self.mpi.world_size == self.mpi.world_rank:
                channel = self.channel[idxij[0]][-1]+self.channel[idxij[1]][-1]
                sta_pair = '{}_{}'.format(f'{self.network[idxij[0]]}_{self.station[idxij[0]]}',
                                          f'{self.network[idxij[1]]}_{self.station[idxij[1]]}')
                # self.create_dataset('COR_{}_{}.h5'.format(sta_pair, channel), idxij, action=action)
                self.logger.info(f'rank {self.mpi.world_rank}: ({i+1}/{len(self.idxij)}) Computing cross-correlation of {self.reftime.strftime('%Y%j%H%M%S')}-{sta_pair}')
                ccf = fftpack.ifft(self.fftarr[idxij[0]]*np.conj(self.fftarr[idxij[1]]), self.nft).real
                ccf = fftpack.ifftshift(ccf)[mid_pos-nlag:mid_pos+nlag+1]
                ff = join('COR_{}_{}.h5'.format(sta_pair, channel))
                self.write_cc(ff, ccf)
        self.mpi.synchronize_all()

    def create_dataset(self, fname, idxij, exist_ok=True, action='w'):
        os.makedirs(self.para.outpath, exist_ok=exist_ok)
        if os.path.isfile(join(self.para.outpath, fname)) and action == 'a':
            return
        with NoiseFile(join(self.para.outpath,fname), action=action) as nf:
            nf.write_header({'network1': self.network[idxij[0]], 'station1': self.station[idxij[0]],
                             'channel1': self.channel[idxij[0]], 'network2': self.network[idxij[1]],
                             'station2': self.station[idxij[1]], 'channel2': self.channel[idxij[1]],
                             'delta': self.dt, 'lag': self.para.maxlag})

    def write_cc(self, fname, data):
        with NoiseFile(join(self.para.outpath,fname), action='a') as nf:
            new_tags = self.reftime.strftime('%Y%m%d%H%M%S')
            nf.write_data(data, name=new_tags)

    def clean(self):
        if self.mpi.world_rank == 0:
            for ff in glob.glob(join(self.para.outpath, '*.h5')):
                os.remove(ff)
        self.mpi.synchronize_all()

    def finalize(self):
        self.mpi.finalize_mpi()
        self.log.destroy_logger()
        self.fft_win.free()
        del self.fftarr
        del self.rawst
        del self.station
        del self.network
        del self.channel


if __name__ == '__main__':
    folder = '/home/xumj/CCCN/example/2016.01.01'
    suf = 'SAC'
    # transf(folder,suf,0.01)
