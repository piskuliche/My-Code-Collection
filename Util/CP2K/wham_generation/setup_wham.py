#!/usr/bin/env python

"""
This is an example code for the calculation of a Potential of Mean Force 
from CP2K or lammps simulations.
This sets up the various windows and gets them ready to run.
Once you run this code - you should submit run_windows.sh
Then you should run Alan Grossfield's WHAM code.
That can be downloaded at http://membrane.urmc.rochester.edu/content/wham
Copyright July 2018 Zeke Piskulich, University of Kansas.

Usage:
    1) Choose # bins and force constants 0.005 for barrier and 0.001 long region is decent
    2) python setup_wham.py 1 to generate wham_metadata.info and run wham file.
    3) python setup_wham.py 2 to run cp2k wham or python setup_wham.py 3 to run lammps wham
    4) sbatch run_windows.sh
    5) bash conv.sh (if cp2k)
    6) Run alan grossfield's wham code as needed
    7) Conversion note: 1 kj/mol/Ang/Ang = 0.0085 hartree/bohr/bohr

Cheers!

"""
import os, shutil, sys
import numpy as np


def generate_forceconsts():
    f = open('force.consts','w')
    num_consts = int(raw_input('How many force consts?\n'))
    prev_bins = 0
    for i in range(num_consts):
        num_bins = int(raw_input('How many bins for const %s?\n' % str(i+1)))
        start = float(raw_input('What is the starting distance (in A)?'))
        end = float(raw_input('What is the ending distance (in A)?'))
        sep = end - start
        kval = float(raw_input('What is the force constant (in hartree/bohr^2 aka cp2k internal)?\n'))
        for j in range(num_bins):
            ro = start+j*sep/float(num_bins)
            f.write('%s %s %s\n' % (prev_bins+1+j, kval, ro))
        prev_bins += num_bins

    f.close()

def setup_cp2k():
    K, ro = np.genfromtxt('force.consts', usecols=(1,2), unpack=True)
    bins = len(K)
    
    angpau=0.529177249
    kcalmolphartree=627.509
    conv = kcalmolphartree/angpau**2
    meta = open('wham_metadata.info', 'w')
    sub = open('run_windows.sh','w')

    sub.write('#!/bin/bash\n')
    sub.write('#SBATCH --job-name=wham_bin\n')
    sub.write('#SBATCH --partition=sixhour\n')
    sub.write('#SBATCH --output=wham_bin-%A_%a.out\n')
    sub.write('#SBATCH --nodes=1\n')
    sub.write('#SBATCH --ntasks-per-node=10\n')
    sub.write('#SBATCH --constraint=intel\n')
    sub.write('#SBATCH --mem=50G\n')
    sub.write('#SBATCH --time=06:00:00\n')
    sub.write('#SBATCH --array 0-%s\n\n\n\n' % (bins))

    sub.write('module load cp2k/6.1/popt\n\n')


    sub.write('cd $SLURM_ARRAY_TASK_ID\n')
    sub.write('mpirun -np 10 cp2k.popt inp_const.cp2k > run1.new\n')
    sub.write("sed 's/\([ \t]\+[^ \t]*\)\{3\}$//' out.colvar > lif.distance\n")
    sub.write('cd ../\n')

    for i in range(0,bins):
        if not os.path.isdir(str(i)):
            os.mkdir(str(i))
        f = open(str(i)+'/collective.inc','w')
        f.write('\t&COLLECTIVE\n')
        f.write('\t COLVAR 1\n')
        f.write('\t INTERMOLECULAR TRUE\n')
        f.write('\t TARGET [angstrom] %s\n' % ro[i])
        f.write('\t\t&RESTRAINT\n')
        f.write('\t\t K %s\n' % K[i])
        f.write('\t\t&END RESTRAINT\n')
        f.write('\t&END COLLECTIVE\n')
        f.close()
        shutil.copyfile('inp_const.cp2k',str(i)+'/inp_const.cp2k')
        #shutil.copyfile('conv.py',str(i)+'/conv.py')
        # Outputs in angstrom, kcal/(mol*angstrom)
        meta.write('%s/output.distance %s %s\n' % (str(i), ro[i] , K[i]*conv))
        
    meta.close()
    sub.close

def setup_lmps():
    K, ro = np.genfromtxt('force.consts', usecols=(1,2), unpack=True)
    bins = len(K)
    angpau=0.529177249
    kcalmolphartree=627.509
    conv = kcalmolphartree/angpau**2
    meta = open('wham_metadata.info', 'w')
    sub = open('run_windows.sh','w')

    sub.write('#SBATCH --job-name=wham_bin\n')
    sub.write('#SBATCH --partition=sixhour\n')
    sub.write('#SBATCH --constraint=intel\n')
    sub.write('#SBATCH --output=output.log\n')
    sub.write('#SBATCH --nodes=1\n')
    sub.write('#SBATCH --ntasks-per-node=10\n')
    sub.write('#SBATCH --time=6:00:00\n')
    sub.write('#SBATCH --array 0-%s\n\n\n' % (bins))

    sub.write('module load lammps/3Mar2020\n\n')


    sub.write('cd $SLURM_ARRAY_TASK_ID\n')
    sub.write('mpirun lmp_mpi <in.ip \n')
    sub.write("sed -e '/TIMESTEP/,+8d' tmp.dump > lif.distance\n")
    sub.write('cd ../\n')

    for i in range(0,bins):
        if not os.path.isdir(str(i)):
            os.mkdir(str(i))
        f = open(str(i)+'/collective.inc','w')
        f.write('bond_coeff  2 %s %s\n' %(K[i]*conv, ro[i]))
        f.close()
        shutil.copyfile('in.ip',str(i)+'/in.ip')
        # Outputs in angstrom, kcal/(mol*angstrom)
        meta.write('%s/lif.distance %s %s\n' % (str(i), ro[i], K[i]*2.0*conv))

    meta.close()
    sub.close

if __name__ == "__main__":
    # runs and chooses the program.
    if len(sys.argv) == 1:
        print("Usage: python setup_wham.py -type [1,2 or 3] ")
        sys.exit("For more information run program with '-h' option")
    if "-h" in sys.argv:
        print("Usage: python setup_wham.py -type [name] ")
        print("This code sets up the wham calculation and sets the force constants for each window.")
        print("type: 1) generate force constants, 2) cp2k setup, 3) lammps setup")
        print("Note run 1 then 2 or 1 then 3 not 1 2 3.")
        sys.exit("Exiting")
    if "-type" in sys.argv:
        index = sys.argv.index("-type")+1
        type = int(sys.argv[index])
    else:
        sys.exit("type must be specified")

    if type == 1:
        generate_forceconsts()
    elif type == 2:
        setup_cp2k()
    elif type == 3:
        setup_lmps()
    else:
        print("Invalid options chosen - try again")

