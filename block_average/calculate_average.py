#!/usr/bin/env python

import numpy as np
import argparse
import math
from argparse import RawTextHelpFormatter,SUPPRESS

parser = argparse.ArgumentParser(description='''Finds average value of a data file''',formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', help ="File Name")
parser.add_argument('-b', help ="Number of Blocks")
parser.add_argument('-s', help ="Number of lines to skip")
parser.add_argument('-n', help ="Name of column")
parser.add_argument('-col', help = "Column number i.e. 0, 1, 2 ..")
args = parser.parse_args()
filename=args.f
colname=args.n
colnum=int(args.col)
block=int(args.b)
skip=int(args.s)
t_val=stats.t.ppf(0.975,block-1)


data = np.genfromtxt(filename,dtype=float,usecols=(colnum),unpack=True, skip_header=skip)
print len(data)
print data
blocksize=len(data)/block
print blocksize
avg = np.zeros(block, dtype = float)
AV=0
for b in range(0,block,1):
    start=b*blocksize
    end=(b+1)*blocksize
    for i in range(start, end, 1):
        avg[b] = avg[b] + data[i]
        print data[i]
    print avg
    avg[b] /= blocksize
    print avg[b]
    AV+=avg[b]
AV= AV / block
if (block > 1):
    STDEV=np.std(avg, ddof=1)*t_val/np.sqrt(block)
else:
    STDEV=0
outfile=colname+".avg"
f=open(outfile,'w')

f.write(str(colname)+" "+str(AV)+" "+str(STDEV)+"\n")
f.close()



