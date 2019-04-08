#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os
"""
This is a python program which creates a box of acn for use in lammps simulations.
"""

if len(sys.argv) != 3:
    print("Usage python build_acn.py numacns type")
    exit()

numacns=int(sys.argv[1])
type=int(sys.argv[2])


def gen_packmol(n):
    L=(n/0.034)**(1/3.)-2
    pmolfile = "acn.pmol"
    pm = open(pmolfile,'w')
    pm.write("tolerance 2.0\n")
    pm.write("filetype xyz\n")
    pm.write("output system.xyz\n")
    pm.write("\n")
    pm.write("structure acn.xyz\n")
    pm.write("  number %s\n" % n)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (L, L, L))
    pm.write("end structure\n")
    pm.close()

    acnfile = "acn.xyz"
    acn = open(acnfile, 'w')
    acn.write("3\n")
    acn.write("acn\n")
    acn.write("ch3          0.000   0.000   0.000\n")
    acn.write("c            1.540   0.000   0.000\n")
    acn.write("n            2.697   0.000   0.000\n")
    acn.close()

    print("Please run the following command: packmol < acn.pmol")


def tolammps(n):
    datafile = "data.acn"
    dat = open(datafile, 'w')
    x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
    q=[-0.269, 0.129, -0.398]
    type=[1,2,3]
    L = (n/0.034)**(1/3.) 
    dat.write("LAMMPS Description\n")
    dat.write("\n")
    dat.write("%s atoms\n" % (n*3))
    dat.write("%s bonds\n" % (n*2))
    dat.write("%s angles\n" % (n*1))
    dat.write("0 dihedrals\n")
    dat.write("0 impropers\n")
    dat.write("\n")
    dat.write("3 atom types\n")
    dat.write("1 bond types\n")
    dat.write("1 angle types\n")
    dat.write("\n")
    dat.write("0.0 %s xlo xhi\n" % L)
    dat.write("0.0 %s  ylo yhi\n" % L)
    dat.write("0.0 %s  zlo zhi\n" % L)
    dat.write("\n")
    dat.write("Masses\n")
    dat.write("\n")
    dat.write("1 15.035\n")
    dat.write("2 12.011\n")
    dat.write("3 14.007\n")
    dat.write("\n")
    dat.write("Atoms\n")
    dat.write("\n")
    for i in range(n):
        for j in range(3):
            aindex=i*3+1+j
            mindex=i+1
            dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, type[j], q[j], x[aindex-1], y[aindex-1], z[aindex-1]))
    dat.write("\n")
    dat.write("Bonds\n")
    dat.write("\n")
    for i in range(n):
        b1index=i*2+1
        b2index=i*2+2
        btype = 1
        oindex =i*3+1
        h1index=i*3+2
        h2index=i*3+3
        dat.write("%s %s %s %s\n" % (b1index,btype,oindex,h1index))
        dat.write("%s %s %s %s\n" % (b2index,btype,oindex,h2index))
    dat.write("\n")
    dat.write("Angles\n")
    dat.write("\n")
    for i in range(n):
        aindex=i+1
        atype=1
        oindex =i*3+1
        h1index=i*3+2
        h2index=i*3+3
        dat.write("%s %s %s %s %s\n" % (aindex, atype, h1index, oindex, h2index))
    dat.close()

if type == 0:
    gen_packmol(numacns)
elif type == 1:
    tolammps(numacns)
elif type == 3:
    gen_packmol(numacns)
    subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < acn.pmol"],shell=True)
    tolammps(numacns)
else:
    print("Error: Incorrect chosen type")
