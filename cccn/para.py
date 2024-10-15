class Para():
    def __init__(self) -> None:
        self.datapath = './'
        self.outpath = './'
        self.freqmin = 0.02
        self.freqmax = 0.2
        self.nsmooth = 10
        self.timeduration = 900
        self.cut_precentatge = 0.01
        self.suffix = 'SAC'
        self.target_dt = None
        self.reftime = 'day'
        self.maxlag = 500

    def bcast(self, mympi):
        mympi.bcast(self.datapath)
        mympi.bcast(self.outpath)
        mympi.bcast(self.freqmin)
        mympi.bcast(self.freqmax)
        mympi.bcast(self.nsmooth)
        mympi.bcast(self.timeduration)
        mympi.bcast(self.cut_precentatge)
        mympi.bcast(self.suffix)
        mympi.bcast(self.target_dt)
        mympi.bcast(self.reftime)
        mympi.bcast(self.maxlag)
        