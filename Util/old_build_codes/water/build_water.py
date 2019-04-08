#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os
"""
This is a python program which creates a box of water for use in lammps simulations.
"""

if len(sys.argv) != 3:
    print("Usage python build_water.py numwaters type")
    exit()

numwaters=int(sys.argv[1])
type=int(sys.argv[2])


def gen_packmol(n):
    L=(n/0.034)**(1/3.)-2
    pmolfile = "water.pmol"
    pm = open(pmolfile,'w')
    pm.write("tolerance 2.0\n")
    pm.write("filetype xyz\n")
    pm.write("output system.xyz\n")
    pm.write("\n")
    pm.write("structure water.xyz\n")
    pm.write("  number %s\n" % n)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (L, L, L))
    pm.write("end structure\n")
    pm.close()

    waterfile = "water.xyz"
    wat = open(waterfile, 'w')
    wat.write("3\n")
    wat.write("water\n")
    wat.write("O           10.203012       7.603800      12.673000\n")
    wat.write("H            9.625597       8.420323      12.673000\n")
    wat.write("H            9.625597       6.787278      12.673000\n")
    wat.close()

    print("Please run the following command: packmol < water.pmol")


def tolammps(n):
    datafile = "data.water"
    dat = open(datafile, 'w')
    x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
    q=[-0.8476, 0.4238, 0.4238]
    type=[1,2,2]
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
    dat.write("1 15.994\n")
    dat.write("2 1.008\n")
    dat.write("3 16.04\n")
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
    gen_packmol(numwaters)
elif type == 1:
    tolammps(numwaters)
elif type == 3:
    gen_packmol(numwaters)
    subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < water.pmol"],shell=True)
    tolammps(numwaters)
else:
    print("Error: Incorrect chosen type")
