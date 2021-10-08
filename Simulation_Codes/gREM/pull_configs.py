#!/usr/bin/env python

import numpy as np
import sys

nlines=int(sys.argv[1])
nfile=int(sys.argv[2])


out=open("step_pull.sh",'w')
with open("log/log.lammps",'r') as f:
    lines=f.readlines()
    line=lines[-1].strip().split()[1:]
    for win in line:
        out.write("tail -n%d %d/dump-%d.dat >> alldumps.lmpstrj\n"%(nlines,int(win),nfile))
out.close()

