import obspy
import os
from os.path import join, exists
import numpy as np
from obspy.signal.util import next_pow_2
from scipy import fftpack
from multiprocessing import Pool as ThreadPool 
from itertools import repeat
from .util import *
from .para import Para
import time
import pyasdf
import glob


class CrossCorrelation():
    def __init__(self, para=None) -> None:
        if para is None:
            self.para = Para()
        else:
            self.para = para
    
    def read_sac(self):
        if exists(self.para.datapath):
            self.rawst = obspy.read(join(self.para.datapath,'*.'+self.para.suffix))
        else:
            self.rawst = obspy.read(self.para.datapath)
        if self.para.target_dt is None:
            self.dt = self.rawst[0].stats.delta
        else:
            self.rawst.resample(1/self.para.target_dt, no_filter=False)
            self.dt = self.para.target_dt

    def perwhiten(self):
        if self.para.reftime == 'day':
            self.reftime = obspy.UTCDateTime(self.rawst[0].stats.starttime.date)
        elif self.para.reftime == 'hour':
            self.reftime = obspy.UTCDateTime(self.rawst[0].stats.starttime.strftime('%Y%m%d%H'))
        self.nft = int(next_pow_2((self.para.timeend - self.para.timestart)/self.dt))
        self.fftst = self.rawst.copy()
        for tr in self.fftst:
            #------- cut waveform ------ 
            cutbtime = self.reftime+self.para.timestart
            cutetime = self.reftime+self.para.timeend
            if tr.stats.starttime > cutbtime or tr.stats.endtime < cutetime:
                continue
            tr.trim(cutbtime, cutetime)
            tr.detrend('linear')
            tr.detrend('constant')
            #----------normalize----------
            if self.para.nsmooth == 0:
                tr.data /= np.abs(tr.data)
            elif self.para.nsmooth > 0:
                tr.data /= smooth(np.abs(tr.data),half_len=self.para.nsmooth)
            else:
                raise ValueError("Half window length must be greater than zero")
            #tr.write(join(folder,"%s.%s.%s.BHZ.norm" % (tr.stats.network, tr.stats.station, tr.stats.location)), "SAC")
            #----------- Whiten -----------
            f1 = 1/(1.5*(1/self.para.freqmin))
            f4 = 1/(0.75*(1/self.para.freqmax))
            (tr.data,tr.stats.delta) = whiten(tr.data, self.nft, self.dt, f1, self.para.freqmin, self.para.freqmax, f4)

    def clean(self):
        for ff in glob.glob(join(self.para.outpath, '*.h5')):
            os.remove(ff)

    def compute_cc(self, idxij, stapair, nts, mid_pos, lag, cor):
        ccf = fftpack.ifft(self.fftst[idxij[0]].data*np.conj(self.fftst[idxij[1]].data), nts).real
        cor.data = fftpack.ifftshift(ccf)[mid_pos-lag:mid_pos+lag+1]
        cor.stats.network = self.fftst[idxij[1]].stats.network
        cor.stats.station = self.fftst[idxij[1]].stats.station
        cor.stats.sac.kevnm = self.fftst[idxij[0]].stats.station
        cor.stats.channel = self.fftst[idxij[0]].stats.channel[-1]+self.fftst[idxij[1]].stats.channel[-1]
        try:
            cor.stats.sac.stla = self.fftst[idxij[1]].stats.sac.stla
            cor.stats.sac.stlo = self.fftst[idxij[1]].stats.sac.stlo
            cor.stats.sac.evla = self.fftst[idxij[0]].stats.sac.stla
            cor.stats.sac.evlo = self.fftst[idxij[0]].stats.sac.stlo
        except:
            pass
        cor.stats.sac.b = -self.para.maxlag
        cor.stats.delta = self.dt
        cor.stats.starttime = self.reftime
        ff = join(self.para.outpath, 'COR_{}_{}.h5'.format(stapair, cor.stats.channel))
        if not os.path.isfile(ff):
            with pyasdf.ASDFDataSet(ff,mpi=False,compression="gzip-3",mode='w') as ds:
                pass
        with pyasdf.ASDFDataSet(ff,mpi=False,compression="gzip-3",mode='a') as ds:
            # try:ds.add_stationxml(inv1) 
            # except Exception: pass 
            new_tags = self.reftime.strftime('%Y%m%d%H%M%S')
            ds.add_waveforms(cor, tag=new_tags)     

    def docc(self):
        pool = ThreadPool(self.para.nnode)
        if not os.path.exists(self.para.outpath):
            os.makedirs(self.para.outpath)
        ns = len(self.fftst)
        nts = self.nft
        nlag = int(self.para.maxlag/self.dt)
        mid_pos = int(nts/2)
    #   tcorr = np.arange(-nts + 1, nts)
    #   dn = np.where(np.abs(tcorr) <= lag)[0]
        cor = self.fftst[0].copy()
        self.sta_pair = []
        idx_lst = []
        for i in np.arange(ns-1):
            for j in np.arange(i+1,ns):
                staname_src = '{}.{}'.format(self.fftst[i].stats.network, self.fftst[i].stats.station)
                staname_sta = '{}.{}'.format(self.fftst[j].stats.network, self.fftst[j].stats.station)
                self.sta_pair.append("{}_{}".format(staname_src, staname_sta))
                idx_lst.append([i, j])
        # t=time.process_time()
        pool.starmap(self.compute_cc, zip(idx_lst, self.sta_pair, repeat(nts), repeat(mid_pos), repeat(nlag), repeat(cor)))
        # print("%d station pair using %d node:" % (len(sta_pair), self.para.nnode), (time.process_time()-t), "s")
        pool.close()
        pool.join()

if __name__ == '__main__':
    folder = '/home/xumj/CCCN/example/2016.01.01'
    suf = 'SAC'
    # transf(folder,suf,0.01)
