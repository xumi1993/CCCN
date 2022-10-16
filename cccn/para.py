class Para():
    def __init__(self) -> None:
        self.datapath = './'
        self.outpath = './'
        self.freqmin = 0.02
        self.freqmax = 0.2
        self.nsmooth = 10
        self.timestart = 800
        self.timeend = 85600
        self.suffix = 'SAC'
        self.target_dt = None
        self.reftime = 'day'
        self.maxlag = 500
        self.nnode = 1