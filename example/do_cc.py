from cccn.noise import CrossCorrelation
from cccn.post_func import PostProcForNoise


def main():
    cc = CrossCorrelation()
    cc.clean()
    for path in ['dataSAC/2008.010', 'dataSAC/2008.011']:
        cc.para.datapath=path
        cc.para.suffix = 'sac'
        cc.para.target_dt = 1
        cc.read_sac()
        cc.perwhiten()
        cc.run_cc()

if __name__ == '__main__':
    main()