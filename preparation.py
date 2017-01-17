#!/usr/bin/env python

import sys
import getopt
import os
import numpy as np
import obspy
from make_one_folder_norm import transf, perwhiten,docc

def Usage():
    print("preparation.py -w<half-length> -f<f1/f2/f3/f4> -d<dt> -c<cut1/cut2> -l<maxlag> -S<suffix> floder_lst")

argv = sys.argv[1:]
for o in argv:
    if os.path.isfile(o):
        folder_lst = o
        break

try:
    opts,args = getopt.getopt(argv, "w:f:S:d:c:l:")
except:
    print('Arguments are not found!')
    sys.exit(1)
if opts == []:
    Usage()
    sys.exit(1)

for op, value in opts:
    if op == "-w":
        wlen = float(value)
    elif op == "-f":
        f1 = float(value.split('/')[0])
        f2 = float(value.split('/')[1])
        f3 = float(value.split('/')[2])
        f4 = float(value.split('/')[3])
    elif op == "-S":
        suffix = value
    elif op == "-d":
        dt = float(value)
    elif op == "-c":
        cuttime1 = float(value.split('/')[0])
        cuttime2 = float(value.split('/')[1])
    elif op == "-l":
        lag = float(value)
    else:
        Usage()
        sys.exit(1)

with open(folder_lst) as flst:
    for folder in flst.readlines():
        folder = folder.strip()
        folder_name = folder.split()[0]
        reftime = obspy.UTCDateTime(folder.split()[1])
        nt = (cuttime2 - cuttime1)/dt
        transf(folder_name, suffix, dt)
        fft_all = perwhiten(folder_name, dt, wlen, cuttime1, cuttime2, reftime, f1,f2,f3,f4)
        if len(fft_all) <= 1:
            print("not enough event in folder %s" % folder_name)
            continue
        docc(folder_name, fft_all,nt,dt,lag, reftime)
