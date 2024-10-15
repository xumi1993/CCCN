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


def create_station_inv(network, station, stla, stlo, channel, stel=0.0, dt=1):
    inv = Inventory(
        # We'll add networks later.
        networks=[],
        # The source should be the id whoever create the file.
        source="{}.{}".format(network, station))
    net = Network(
        # This is the network code according to the SEED standard.
        code=network,
        # A list of stations. We'll add one later.
        stations=[])
    sta = Station(
        code=station,
        latitude=stla,
        longitude=stlo,
        elevation=stel)
    cha = Channel(
        # This is the channel code according to the SEED standard.
        code=channel,
        # This is the location code according to the SEED standard.
        location_code="",
        # Note that these coordinates can differ from the station coordinates.
        latitude=stla,
        longitude=stlo,
        elevation=stel,
        depth=0.0,
        azimuth=0.0,
        dip=-90.0,
        sample_rate=1/dt)
    sta.channels.append(cha)
    net.stations.append(sta)
    inv.networks.append(net)
    return inv


def create_quake(network, station, evla, evlo, evel=0.0):
    ori = Origin(longitude=evlo,
                latitude=evla,
                depth=-evel,
                # comments='{}.{}'.format(network, station),
                )
    evt = Event(resource_id=ResourceIdentifier('{}.{}'.format(network, station)),
                origins=[])
    evt.origins.append(ori)
    return evt            


class CrossCorrelation():
    def __init__(self, para=None) -> None:
        self.mpi = MyMPI()
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
    
    def read_sac(self):
        if self.mpi.world_rank == 0:
            if exists(self.para.datapath):
                fname = join(self.para.datapath,'*.'+self.para.suffix)
            else:
                fname = self.para.datapath
            try:
                self.rawst = obspy.read(fname)
            except:
                return False
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
        self.staloc, self.sta_win = self.mpi.prepare_shm([self.ntr, 3], np.float64)
        if self.mpi.world_rank == 0:
            for i, tr in enumerate(self.rawst):
                stel = tr.stats.sac.stel if hasattr( tr.stats.sac, 'stel') else 0.0
                self.network.append(tr.stats.network)
                self.station.append(tr.stats.station)
                self.channel.append(tr.stats.channel)
                self.staloc[i] = [tr.stats.sac.stla, tr.stats.sac.stlo, stel]
        self.network = self.mpi.bcast(self.network)
        self.station = self.mpi.bcast(self.station)
        self.channel = self.mpi.bcast(self.channel)
        self.mpi.sync_from_main(self.staloc)


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
            for tr in fftst:
                #------- cut waveform ------ 
                cutbtime = self.reftime+timestart
                cutetime = self.reftime+timeend
                if tr.stats.starttime > cutbtime or tr.stats.endtime < cutetime:
                    fftst.remove(tr)
                    self.rawst.remove(tr)
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

    def prepare_cc(self):
        self.idxij = []
        for i in np.arange(self.ntr-1):
            for j in np.arange(i+1, self.ntr):
                self.idxij.append([i, j])
        self.mpi.synchronize_all()

    def run_cc(self):
        self.prepare_cc()
        mid_pos = int(self.nft/2)
        nlag = int(self.para.maxlag/self.dt)
        for i, idxij in enumerate(self.idxij):
            if i % self.mpi.world_size == self.mpi.world_rank:
                ccf = fftpack.ifft(self.fftarr[idxij[0]]*np.conj(self.fftarr[idxij[1]]), self.nft).real
                ccf = fftpack.ifftshift(ccf)[mid_pos-nlag:mid_pos+nlag+1]
                channel = self.channel[idxij[0]][-1]+self.channel[idxij[1]][-1]
                cor = SACTrace()
                cor.data = ccf
                cor.kcmpnm = channel
                cor.kevnm = f'{self.network[idxij[0]]}.{self.station[idxij[0]]}'
                try:
                    cor.stla = self.staloc[idxij[1]][0]
                    cor.stlo = self.staloc[idxij[1]][1]
                    cor.evla = self.staloc[idxij[0]][0]
                    cor.evlo = self.staloc[idxij[0]][1]
                except:
                    pass
                cor.b = -self.para.maxlag
                cor.delta = self.dt
                cor_tr = cor.to_obspy_trace()
                cor_tr.stats.starttime = self.reftime
                sta_pair = '{}_{}'.format(f'{self.network[idxij[0]]}.{self.station[idxij[0]]}',
                                          f'{self.network[idxij[1]]}.{self.station[idxij[1]]}')
                ff = join('COR_{}_{}.h5'.format(sta_pair, channel))
                self.create_dataset(ff, idxij)
                self.write_cc(ff, cor_tr)

    def create_dataset(self, fname, idxij, exist_ok=True):
        os.makedirs(self.para.outpath, exist_ok=exist_ok)
        if os.path.isfile(join(self.para.outpath, fname)):
            return
        with pyasdf.ASDFDataSet(join(self.para.outpath,fname), mpi=False,compression="gzip-3",mode='w') as ds:
            stel = self.staloc[idxij[1]][2]
            evel = self.staloc[idxij[0]][2]
            inv = create_station_inv(self.network[idxij[1]], self.station[idxij[1]],
                                        self.staloc[idxij[1]][0],
                                        self.staloc[idxij[1]][1],
                                        self.channel[idxij[1]],
                                        stel,
                                        dt=self.dt)
            quake = create_quake(self.network[idxij[0]],
                                 self.station[idxij[0]],
                                 self.staloc[idxij[0]][0],
                                 self.staloc[idxij[0]][1],
                                 evel)
            ds.add_stationxml(inv)
            ds.add_quakeml(quake)

    def write_cc(self, fname, tr):
        with pyasdf.ASDFDataSet(fname, mpi=False,compression="gzip-3",mode='a') as ds:
            new_tags = self.reftime.strftime('%Y%m%d%H%M%S')
            ds.add_waveforms(tr, tag=new_tags)  
    

    def clean(self):
        if self.mpi.world_rank == 0:
            for ff in glob.glob(join(self.para.outpath, '*.h5')):
                os.remove(ff)

    # def compute_cc(self, idxij, stapair, nts, mid_pos, lag, cor):
    #     ccf = fftpack.ifft(self.fftst[idxij[0]].data*np.conj(self.fftst[idxij[1]].data), nts).real
    #     cor.data = fftpack.ifftshift(ccf)[mid_pos-lag:mid_pos+lag+1]
    #     cor.stats.network = self.fftst[idxij[1]].stats.network
    #     cor.stats.station = self.fftst[idxij[1]].stats.station
    #     cor.stats.sac.kevnm = self.fftst[idxij[0]].stats.station
    #     cor.stats.channel = self.fftst[idxij[0]].stats.channel[-1]+self.fftst[idxij[1]].stats.channel[-1]
    #     try:
    #         cor.stats.sac.stla = self.fftst[idxij[1]].stats.sac.stla
    #         cor.stats.sac.stlo = self.fftst[idxij[1]].stats.sac.stlo
    #         cor.stats.sac.evla = self.fftst[idxij[0]].stats.sac.stla
    #         cor.stats.sac.evlo = self.fftst[idxij[0]].stats.sac.stlo
    #     except:
    #         pass
    #     cor.stats.sac.b = -self.para.maxlag
    #     cor.stats.delta = self.dt
    #     cor.stats.starttime = self.reftime
    #     ff = join(self.para.outpath, 'COR_{}_{}.h5'.format(stapair, cor.stats.channel))
    #     if not os.path.isfile(ff):
    #         with pyasdf.ASDFDataSet(ff,mpi=False,compression="gzip-3",mode='w') as ds:
    #             if hasattr( self.fftst[idxij[1]].stats.sac, 'stel'):
    #                 stel = self.fftst[idxij[1]].stats.sac.stel
    #             else:
    #                 stel = 0.0
    #             if hasattr(self.fftst[idxij[0]].stats.sac, 'stel'):
    #                 evel = self.fftst[idxij[0]].stats.sac.stel
    #             else:
    #                 evel = 0.0
    #             inv = create_station_inv(cor.stats.network, cor.stats.station,
    #                                      self.fftst[idxij[1]].stats.sac.stla,
    #                                      self.fftst[idxij[1]].stats.sac.stlo,
    #                                      cor.stats.channel,
    #                                      stel,
    #                                      dt=self.fftst[idxij[1]].stats.delta)
    #             quake = create_quake(self.fftst[idxij[0]].stats.network,
    #                                  self.fftst[idxij[0]].stats.station,
    #                                  self.fftst[idxij[0]].stats.sac.stla,
    #                                  self.fftst[idxij[0]].stats.sac.stlo,
    #                                  evel)
    #             ds.add_stationxml(inv)
    #             ds.add_quakeml(quake)
    #     with pyasdf.ASDFDataSet(ff,mpi=False,compression="gzip-3",mode='a') as ds:
    #         # try:ds.add_stationxml(inv1) 
    #         # except Exception: pass 
    #         new_tags = self.reftime.strftime('%Y%m%d%H%M%S')
    #         ds.add_waveforms(cor, tag=new_tags)     



if __name__ == '__main__':
    folder = '/home/xumj/CCCN/example/2016.01.01'
    suf = 'SAC'
    # transf(folder,suf,0.01)
