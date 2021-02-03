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
    tmp, mols = zip(*sorted(zip(id,mols)))
    tmp, typ = zip(*sorted(zip(id,typ)))
    id, Q   = zip(*sorted(zip(id,Q)))
    np.savetxt('molinfo.dat', np.c_[id,mols,typ,Q], fmt=('%d','%d','%d', '%1.4f'))
    print('Molecular information saved to molinfo.dat')
    print('There are %d atoms and %d items in molinfo.dat' % (natoms, len(mols)))
    return natoms, L

def read_log(filename, ecol, volcol):
    with open(filename, 'r') as f:
        lines=f.readlines()
        flag = 0
        for line in lines:
            if "Loop" in line.split():
                flag = 0
            if flag == 1:
                e.append(float(line.split()[ecol]))
                vol.append(float(line.split()[volcol]))
            if "Step" in line.split():
                flag =1
                e = []
                vol = []
        Lval = np.asarray(vol)**(1/3.)
        np.savetxt('e_init.out', np.c_[e])
        np.savetxt('vol_init.out', np.c_[vol])
        np.savetxt('L.dat', np.c_[Lval])

    
if __name__ == "__main__":
    import sys,argparse
    realL=-1

    # Read in input parameters
    parser=argparse.ArgumentParser()
    parser.add_argument('-f', default="data.lmps", help='Data file name')
    parser.add_argument('-dr', default=0.1, help='Shell thickness')
    parser.add_argument('-s1', default=1, help='First atom selection')
    parser.add_argument('-s2', default=1, help='Second atom selection')
    parser.add_argument('-sconfig', default=0, help='Starting frame for calculation')
    parser.add_argument('-econfig', default=1000, help='Ending frame for calculation')
    parser.add_argument('-or_int', default=1, help='Every nth frame should be used')
    parser.add_argument('-sskip', default=2, help='Number of lines to skip at frame start')
    parser.add_argument('-eskip', default=0, help='Number of lines to skip at frame end')
    parser.add_argument('-realL', default=-5.0, help='Box length in Angstroms')
    parser.add_argument('-ecol', default=2, help='Location of energy in log file')
    parser.add_argument('-volcol', default=10, help='Location of volume in log file')
    parser.add_argument('-logname', default="log.lammps", help='Name of log file')
    args = parser.parse_args()

    filename=str(args.f)
    dr=float(args.dr)
    selec1=int(args.s1)
    selec2=int(args.s2)
    startconfig=int(args.sconfig)
    endconfig=int(args.econfig)
    or_int=int(args.or_int)
    startskip=int(args.sskip)
    endskip=int(args.eskip)
    realL=float(args.realL)
    ecol=int(args.ecol)
    volcol=int(args.volcol)
    logname=str(args.logname)
    


    natoms, L = read_mols(filename)
    if realL > 0:
        L = realL
    create_template(natoms,L,selec1,selec2, startconfig, endconfig, or_int, dr,startskip, endskip)
    read_log(logname, ecol, volcol)
