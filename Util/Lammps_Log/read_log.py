#!/usr/bin/env python


import numpy as np
import argparse
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument('-f', default="log.lammps", help='Log file to parse')
parser.add_argument('-nb', default=5, help='Number of blocks')
parser.add_argument('-info', default="value", help='Information to append to output files')


args = parser.parse_args()
fname = str(args.f)
nblocks=int(args.nb)
info=str(args.info)

readflag=0
runs = []
runstep = -1
items = []
with open(fname, 'r') as f:
    lines=f.readlines()
    for line in lines:
        if "Step" in line:
            readflag = 1
            runstep += 1
            data = {}
            if runstep == 0:
                for item in line.split():
                    items.append(item)
            for item in items:
                data[item] = []
        if "Loop time" in line:
            readflag = 0
            runs.append(data)
        if readflag == 1:
            if line.split()[0] != items[0]:
                for i in range(len(items)):
                    data[items[i]].append(float(line.split()[i]))

print("There are %s run commands" % len(runs))
for r in range(len(runs)):
    print("*****************************************************************")
    print("Run %d has %d data points" % (r, len(runs[r][items[0]])))
    if runs[r]["Volume"][0]==runs[r]["Volume"][1]:
        print("Run Type is NVT")
    else:
        print("Run Type is NPT")

    nitems = len(runs[r][items[0]])-1
    nperb = int(nitems/nblocks)
    t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

    for item in items:
        if r == len(runs)-1:
            np.savetxt(str(item)+'_init.out', np.c_[runs[r][item]])
        bl = []
        for b in range(nblocks):
            bstart = b*nperb + 1
            bend = bstart + nperb
            bl.append(np.average(runs[r][item][bstart:bend]))
        std = np.std(bl)*t_val
            
        print("The average of item %s is %.5f +/- %.5f" %(item,np.average(runs[r][item]),std))
        f = open("run"+str(r+1)+"_"+str(item)+'.avg','w')
        f.write('%s %.5f %.5f\n' % (info,np.average(runs[r][item]),std))
        f.close()
