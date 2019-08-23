#!/usr/bin/env python
"""
This code does a few things: 
    1) it reads the data file to create a mols.dat
    2) it creates a blank gofr input template
"""
import numpy as np

def create_template(natoms,L,selec1,selec2, startconfig, endconfig, or_int, dr,sskip, eskip):
    f = open('gofr.inp', 'w')
    print("Creating gofr.inp")
    f.write('&nml\n')
    f.write('nfile="traj.xyz",\n')
    f.write('molfile="molinfo.dat",\n')
    f.write('L=%1.6f,\n' % np.abs(L))
    f.write('natoms=%d,\n' % natoms)
    f.write('startconfig=%d,\n' % startconfig)
    f.write('endconfig=%d,\n' % endconfig)
    f.write('or_int=%d,\n' % or_int)
    f.write('startskip=%d,\n'% sskip)
    f.write('endskip=%d,\n'% eskip)
    f.write('nblocks=5,\n')
    f.write('selec1=%d,\n' % selec1)
    f.write('selec2=%d,\n' % selec2)
    f.write('dr=%1.5f\n' % dr)
    f.write('/\n')
    f.close()
    return

def read_mols(filename):
    savelines=[]
    counts = 0
    sskip = 100
    eskip = 0
    elinesf = 0
    natoms=0
    with open(filename, 'r') as f:
        lines=f.readlines()
        for line in lines:
            counts += 1
            if "atoms" in line:
                natoms=int(line.split()[0])
            # Sets the number of lines for skipping
            if "xlo" in line:
                L = float(line.split()[0])-float(line.split()[1])
            if "Atoms" in line:
                sskip = counts
                eskip = len(lines)-counts-natoms-1
            # Have to remove empty lines from the footer because
            # genfromtxt is dumb and ignores empty lines
            if counts > sskip + natoms:
                if line.strip()=="":
                    eskip-=1
    id,mols,typ,Q = np.genfromtxt(filename, usecols=(0,1,2,3), unpack=True, skip_header=sskip, skip_footer=eskip)
    np.savetxt('molinfo.dat', np.c_[id,mols,typ,Q], fmt=('%d','%d','%d', '%1.4f'))
    print('Molecular information saved to molinfo.dat')
    print('There are %d atoms and %d items in molinfo.dat' % (natoms, len(mols)))
    return natoms, L

if __name__ == "__main__":
    import sys
    realL=-1
    filename = "-h"
    if len(sys.argv) < 3:
        print("Usage: setup_gofr.py filename dr [selec1 selec2] [startconfig endconfig or_int] [startskip endskip]")
        sys.exit()
    filename = str(sys.argv[1])
    if filename == "-h":
        print("Usage: setup_gofr.py filename dr [selec1 selec2] [startconfig endconfig or_int] [startskip endskip]")
        sys.exit()
    dr = float(sys.argv[2])
    if len(sys.argv) > 3 and len(sys.argv) <= 5:
        selec1 = int(sys.argv[3])
        selec2 = int(sys.argv[4])
        startconfig=0
        endconfig=1000
        or_int=1
        startskip=2
        endskip=0
    elif len(sys.argv) > 5 and len(sys.argv) <= 8:
        selec1 = int(sys.argv[3])
        selec2 = int(sys.argv[4])
        startconfig = int(sys.argv[5])
        endconfig = int(sys.argv[6])
        or_int = int(sys.argv[7])
        startskip=2
        endskip=0
    elif len(sys.argv) > 8:
        selec1 = int(sys.argv[3])
        selec2 = int(sys.argv[4])
        startconfig = int(sys.argv[5])
        endconfig = int(sys.argv[6])
        or_int = int(sys.argv[7])
        startskip = int(sys.argv[8])
        endskip = int(sys.argv[9])
        realL = float(sys.argv[10])
    else: 
        selec1 = 1
        selec2 = 1
        startconfig=0
        endconfig=1000
        or_int=1
        startskip=2
        endskip=0
    natoms, L = read_mols(filename)
    if realL > 0:
        L = realL
    create_template(natoms,L,selec1,selec2, startconfig, endconfig, or_int, dr,startskip, endskip)
