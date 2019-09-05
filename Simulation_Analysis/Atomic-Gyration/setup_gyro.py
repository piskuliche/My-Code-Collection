#!/usr/bin/env python
"""
This code does a few things: 
    1) it reads the data file to create a mols.dat
    2) it creates a blank Rg input template
"""
import numpy as np

def create_template(natoms,L,selec, startconfig, endconfig, or_int, sskip, eskip,nblocks):
    from scipy import stats
    t_val=stats.t.ppf(0.975,nblocks-1)
    f = open('gyro.inp', 'w')
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
    f.write('tval=%s,\n' % t_val)
    f.write('nblocks=%d,\n' % nblocks)
    f.write('selec=%d,\n' % selec)
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
    filename = "-h"
    if len(sys.argv) < 3:
        print("Usage: setup_gyro.py filename [selec] [startconfig endconfig or_int] [startskip endskip nblocks]")
        sys.exit()
    filename = str(sys.argv[1])
    realL=-1.0
    if filename == "-h":
        print("Usage: setup_gyro.py filename [selec] [startconfig endconfig or_int] [startskip endskip nblocks]")
        sys.exit()
    if len(sys.argv) > 2 and len(sys.argv) <= 3:
        selec = int(sys.argv[2])
        startconfig=0
        endconfig=1000
        or_int=1
        startskip=2
        endskip=0
        nblcoks=5
    elif len(sys.argv) > 3 and len(sys.argv) <= 6:
        selec = int(sys.argv[2])
        startconfig = int(sys.argv[3])
        endconfig = int(sys.argv[4])
        or_int = int(sys.argv[5])
        startskip=2
        endskip=0
        nblocks=5
    elif len(sys.argv) > 6:
        selec = int(sys.argv[2])
        startconfig = int(sys.argv[3])
        endconfig = int(sys.argv[4])
        or_int = int(sys.argv[5])
        startskip = int(sys.argv[6])
        endskip = int(sys.argv[7])
        nblocks = int(sys.argv[8])
        realL = float(sys.argv[9])
    else: 
        selec = 1
        startconfig=0
        endconfig=1000
        or_int=1
        startskip=2
        endskip=0
        nblocks=5
    
    natoms, L = read_mols(filename)
    if realL > 0:
        L=realL
    create_template(natoms,L,selec, startconfig, endconfig, or_int, startskip, endskip,nblocks)
