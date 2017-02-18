#!/usr/bin/env python

import sys
import getopt
import os
import numpy as np
import obspy
import time
from cccn.make_one_folder import transf, transf_hinet, perwhiten,docc

def Usage():
    print("preparation.py -f<f1/f2/f3/f4> -d<dt> -c<cut1/cut2> -l<maxlag>")
    print("[-n<node>] [-C<channel_list>] [-t] [-w<half-length>] [-S<suffix>]")
    print("[-F<seed|hinet>] floder_lst")

global fft_all
argv = sys.argv[1:]
for o in argv:
    if os.path.isfile(o):
        folder_lst = o
        break

try:
    opts,args = getopt.getopt(argv, "w:f:S:d:c:l:t:C:n:F:")
except:
    print('Arguments are not found!')
    sys.exit(1)
if opts == []:
    Usage()
    sys.exit(1)

suffix = "SAC"
wlen = 0
istransf = True
ch = ['Z']
node = 1
dformat = "seed"
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
    elif op == "-t":
        istransf = False
    elif op == "-C":
        ch = [cha for cha in value]
    elif op == "-n":
        node = int(value)
    elif op == "-F":
        if value in ["seed", "hinet"]:
            dformat = value
        else:
            print("wrong option in \"-F\" argument")
            sys.exit(1)
    else:
        Usage()
        sys.exit(1)

with open(folder_lst) as flst:
    for folder in flst.readlines():
        folder = folder.strip()
        print(folder)
        folder_name = folder.split()[0]
        reftime = obspy.UTCDateTime(folder.split()[1])
        nt = int(np.floor((cuttime2 - cuttime1)/dt))
        if istransf:
            if dformat == "seed":
                transf(folder_name, suffix, dt,ch=ch)
            elif dformat == "hinet":
                transf_hinet(folder_name, suffix, dt,ch=ch)
            else:
                print("wrong option in \"-F\" argument")
                sys.exit(1)
        ev_num = perwhiten(folder_name, dt, wlen, cuttime1, cuttime2, reftime, f1,f2,f3,f4,ch=ch)
        if ev_num <= 1:
            print("not enough event in folder %s" % folder_name)
            continue
        sta_pair = docc(folder_name ,nt,dt,lag, reftime,f2,f3, node)
