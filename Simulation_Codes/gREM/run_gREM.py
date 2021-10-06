#!/usr/bin/env python
########################################################################################
# This code is a personal version written by Zeke Piskulich
# It takes logfiles from replica exchange simulations and then pulls the data from them
# it takes heavy inspiration from the code written by David Stelter for the gREM method
# in 2015. zzz
# For questions/concerns reach out to Zeke Piskulich (piskulichz@gmail.com)
# Copyright 2021 - Boston University
########################################################################################

import numpy as np
import argparse
import os,pickle,sys
from walker import walkers
from stwham import run_stwham
from scipy import stats
from sortdumps import frame, get_frames,find_leaflet

consts={"kb":0.0019872041} # kcal/mol

# Read Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-setup', default=0, type=int, help='Should it set up variables [default=0 (off)]')
parser.add_argument('-start', default=0,type=int, help='Starting index for calculating [default=0]')
parser.add_argument('-stop', default=5, type=int, help='Ending index for calculating [default=5]')
parser.add_argument('-workdir', default=".", type=str, help='Working directory location [default=./]')
parser.add_argument('-eta', default=-0.03, type=float, help='Eta value used in the simulations [default=-0.03]')
parser.add_argument('-Ho', default=-30000, type=float, help='H0 value used in the simulations [default=-30000]')
parser.add_argument('-nbins', default=20, type=int, help='Number of bins used for histogramming [default=20]')
parser.add_argument('-nb', default=5, type=int, help='Number of blocks used for block averaging [default=5]')
parser.add_argument('-walkdown', default=0, type=int, help='Runs walkdown instead of analysis [default=0 (off)]')
parser.add_argument('-nprocs', default=4, type=int, help='(Walkdown Only) Number of processors [default=4]')
parser.add_argument('-lmp',default="/share/pkg.7/lammps/3Mar2020/install/bin/lmp_mpi", type=str, help='(Walkdown Only) Path to LAMMPS Executable [default=/share/pkg.7/lammps/3Mar2020/install/bin/lmp_mpi]')
parser.add_argument('-inp', default="start.gREM", type=str, help='Lammps input file name [default=start.gREM]')
parser.add_argument('-lambdafile', default="grem.include", type=str, help='location of lambda information [default=grem.include]')
parser.add_argument('-restart', default=0, type=int, help='Should it read a restart file instead of going through reading files? [default=0]')
parser.add_argument('-dumpbase', default="None", type=str, help='Base name of the dump file i.e. dump in dump-1.dat. [default=None]')
parser.add_argument('-safe', default=1, type=int, help='Determines whether to identify leaflets every step or not [default=1 : every]')
parser.add_argument('-hatom', default=4, type=int, help='Atom type for header atom')
parser.add_argument('-rcut', default=12, type=float, help='Cutoff distance in A')
parser.add_argument('-slurm', default=0, type=int, help='[0] not slurm [1] slurm')
parser.add_argument('-tag', default="dppc", type=str, help='Tag output files [default=dppc]')
args= parser.parse_args()

setup   = args.setup
fstart  = args.start
fend    = args.stop
workdir = args.workdir
eta     = args.eta
Ho      = args.Ho
nbins   = args.nbins
dowalkdown = args.walkdown
nprocs  = args.nprocs
lmp     = args.lmp
lmpin   = args.inp
lambdafile = args.lambdafile
shouldrestart=args.restart
dumpbase = args.dumpbase
safeleaf = args.safe
nblocks  = args.nb
hatom   = args.hatom
rcut    = args.rcut
slurm   = args.slurm
tag     = args.tag
t_value = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

if shouldrestart == 0:
    print("Restart file will be generated!")
elif shouldrestart == 1:
    print("Reading restart file!")
else:
    exit("Invalid option for restart, please select 0 or 1")

if setup == 1: # This does setup for the simulation, generates folders, etc
    print("To Setup run, use start to choose lowest lambda, stop to choose highest, and nbins to choose number of windows")
    # Evenly choose lambdas and bins
    lambdas=np.linspace(fstart,fend,nbins)
    reps=np.linspace(0,nbins-1,nbins)
    # write lammps variables to a file
    f=open(workdir+"/grem.include",'w')
    f.write("variable lambda world ")
    for i in range(nbins):
        f.write("%d " %lambdas[i])
    f.write("\n")
    f.write("variable walker world ")
    # Makes subdirs and write reps to lammps var file
    for i in range(nbins):
        f.write("%d " % reps[i])
        try:
            os.makedirs(str(i))
        except:
            continue
    f.write("\n")
    f.close()
    f=open(workdir+"/last",'w')
    f.write("0\n")
    f.close()
    try:
        os.mkdir(workdir+"/log")
    except:
        print("log folder exists")
    # Setup step has finished, exits.
    exit("Setup complete")


class Logger(object):
    # Note - this class is for writing code output to both screen and a logfile (in case you want to get fancy)
    # Note - this code was taken from 
    # https://stackoverflow.com/questions/14906764/how-to-redirect-stdout-to-both-file-and-console-with-scripting
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("grem-logfile.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    


       
def read_walker(walklog,start,stop):
    """
    Reads data from the walker log file, saving the data in a numpy array
    """
    walkdata=np.genfromtxt(walklog+"-"+str(start), skip_header=3, skip_footer=1,dtype=int)
    for run in range(start+1,stop):
        walkdata=np.append(walkdata,np.genfromtxt(walklog+"-"+str(run),skip_header=3,skip_footer=1,dtype=int),axis=0)
    return walkdata

def read_lambdas(lammpsfile):
    """
    Function that reads the lambdas from a file and saves them.
    """
    vals=[]

    with open(lammpsfile,'r') as f:
        lines=f.readlines()
        for line in lines:
            if "variable lambda world" in line:
                vals = np.array(line.strip().split()[3:],dtype='float')
    np.savetxt("lambdas.dat", np.c_[vals])
    return vals

def gen_walkdown(lambdas):
    """
    This function generates a submit script for the SCC that runs the walkdown when you submit it.

    Note that this might need to be modified based off of cluster environment.
    """
    f = open("walkdown.sh",'w')
    f.write("#!/bin/bash\n")
    f.write("#\n")
    f.write("#$ -l h_rt=11:45:00\n")
    f.write("#$ -j y\n")
    f.write("#$ -N %s-wd\n"%tag)
    f.write("#$ -pe mpi_%d_tasks_per_node %d\n"% (nprocs,nprocs))
    f.write("#$ -V\n")
    f.write("\n")
    f.write("module load gcc/8.3.0 \nmodule load openmpi/3.1.4_gnu-8.1 \nmodule load fftw/3.3.8\nmodule load lammps/3Mar2020\nmodule load codecol/grem\n")

    f.write("NSLOTS=%d\n" % nprocs)

    f.write("run_gREM.py -nprocs %d -lmp %s -inp %s -walkdown 2 -workdir %s -eta %s -Ho %s\n" % (nprocs, lmp, lmpin,workdir,eta,Ho))
    f.write("exit 0\n")
    f.close()
    print("walkdown.sh created in present directory")
    print("to continue: type qsub walkdown.sh")


def walkdown(lambdas):
    """
    This function performs the walkdown using OS system calls
    It actually is called within the submit script created in gen_walkdown
    and then it runs the lammps trajectories moving from walker 0 to walker M in order
    """
    import os
    nreps=len(lambdas)
    print("There are %d replicas" % nreps)
    os.chdir(workdir)
    for rep in range(nreps):
        print("walking replica %d" % rep)
        os.chdir("./%s" % rep)
        if slurm == 0:
            os.system("mpirun -np %d " % (nprocs/2) + lmp + " -sf omp -pk omp 2 -in " + "../" + lmpin + " -var lambda %g -var eta %g -var H0 %g > output" % (lambdas[rep], eta, Ho))
        else:
            os.system("mpirun -np %d " % (nprocs/2) + lmp + " -in " + "../" + lmpin + " -var lambda %g -var eta %g -var H0 %g > output" % (lambdas[rep], eta, Ho))
        os.system("mv final_restart_file final_restart_file0")
        os.system("cp final_restart_file0 restart_file")
        os.system("cp final_restart_file0 ../%s/restart_file" % (rep+1))
        os.chdir("../") # returns to original folder

def calculate_acceptance(walkloc):
    # This calculates how many exchanges are accepted for each replica.
    print(np.shape(walkloc))
    nswp = np.shape(walkloc)[0]-1
    nrep = np.shape(walkloc)[1]
    c=0
    told=[]
    total=[]
    vals={}
    for t in walkloc:
        if c == 0:
            told=t
            for i in told:
                vals[i]=0
            c=c+1
        else:
            msk=(t!=told)*1
            b=0
            for i in t:
                vals[i]=vals[i]+msk[b]
                b=b+1
            h=np.sum(msk)
            total.append(h/nrep)
            told=t
            c = c + 1
    total = np.sum(total)/nswp
    for i in vals:
        print(i, vals[i]/nswp)
    print("total",total)


if __name__ == "__main__":
    sys.stdout = Logger()
    lambdas=read_lambdas(workdir+"/"+lambdafile)
    if dowalkdown == 1:
        gen_walkdown(lambdas)
    elif dowalkdown == 2:
        walkdown(lambdas)
    elif dowalkdown == 0 and shouldrestart == 0:
        # read in the walker data (which window it is in)
        walkloc=read_walker(workdir+"/log/log.lammps",fstart,fend)[:,1:]
        calculate_acceptance(walkloc)
        # read input file
        allwalkers=walkers(workdir,fstart,fend,walkloc,eta,Ho)
        pickle.dump(walkloc,open(workdir+'/walkloc.pckl','wb'))
        pickle.dump(allwalkers,open(workdir+'/allwalkers.pckl','wb'))
    elif dowalkdown == 0 and shouldrestart == 1:
        allwalkers=pickle.load(open(workdir+'/allwalkers.pckl','rb'),encoding='latin1')
        allwalkers.post_process(nbins,nblocks)
        # Block Averaging
        bl_T, bl_S =[], []
        for blk in range(nblocks):
            T_tmp, S_tmp, _, _, _= run_stwham(allwalkers.bl_hist["PotEng"][blk],allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"],allwalkers.lambdas,allwalkers.eta,allwalkers.Ho,blk)
            bl_T.append(T_tmp)
            bl_S.append(S_tmp)
        T_err = np.std(bl_T,axis=0)*t_value
        S_err = np.std(bl_S,axis=0)*t_value
        Th, Sh, bs, be, binsize = run_stwham(allwalkers.hist["PotEng"],allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"],allwalkers.lambdas,allwalkers.eta,allwalkers.Ho,-1)
        for blk in range(nblocks):
            with open("TandS_bl_"+str(blk)+".out",'w') as bout:
                for i in range(bs,be):
                    bout.write("%f %f %f \n" % (allwalkers.minval["PotEng"]+(i*binsize),bl_T[blk][i],bl_S[blk][i]))
        with open("TandS_STWHAM.out",'w') as tout:
            for i in range(bs,be):
                for blk in range(nblocks):
                    if bl_T[blk][i]==0:
                        T_err[i],S_err[i]=0,0
                tout.write("%f %f %f %f %f\n" % (allwalkers.minval["PotEng"]+(i*binsize), Th[i], T_err[i],Sh[i],S_err[i]))
