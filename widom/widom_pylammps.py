import numpy as np
import sys,random,math
from mpi4py import MPI

"""
         bolz : e^-du/kbT array over ninsert
      Bi_step : step value of V*e^-du/kbT
           Bi : sum over steps of V*e^-du/kbT
    MUex_step : step value of e^-du/kbT
         MUex : Excess chemical potential 
       E_step : The energy of the step (before insertions)
"""

def run(steps):
    lmp.command("run %s" % int(steps))

# UNITS BLOCK
kb=0.0019872041 #kcal/molK
T=298.15 # K
beta = (1.0/(kb*T))

# MPI SECTION
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Argument Parser
if len (sys.argv) != 4:
    print("Usage: widom_lammps.py input ninsert nsteps")
    sys.exit()
infile=sys.argv[1]
ninsert = int(sys.argv[2])
nsteps = int(sys.argv[3])

# Random Number Seed
random.seed(27413)

def randnum():
    rand = val = ''.join(random.sample("0123456789", 6))
    return rand

# Import LAMMPS and initialize it
from lammps import PyLammps
lmp = PyLammps()

# read in input file
lmp.file(infile)

#lmp.variable("e equal pe")
run(100)
# calculate the energy of the system
if rank == 0:
    natoms = lmp.system.natoms
    inite = lmp.eval("pe") / natoms
    vol = lmp.system.xhi**3
    print("%s atoms in the system" % natoms)
    print("%s is the initial energy (kcal/mol)" % inite)

lmp.molecule("methane ch4.txt toff 2")
run(100)

Bi_step, Bi, MUex_step, MUex, MUinf = 0.0, 0.0, 0.0, 0.0, 0.0
av_vol, av_rho = 0.0, 0.0
E_step = 0.0

# Open logfile
logfile = "log.widom"
log = open(logfile, 'w')
if rank == 0:
    log.write("Step MUinf MUex\n")
# Begin step loop
for i in range(nsteps):
    print i
    # Zeros the energy values.
    bolz = []
    Bi_step, MUex_step = 0.0, 0.0
    # Runs 1000 steps, computes pe, grabs volume, grabs density
    run(1000)
    if rank == 0:
        E_step=lmp.eval("pe")
        vol = lmp.system.xhi**3
        rho = lmp.eval("density")
    # Loop over insertion attempts.
    for j in range(ninsert):
        print j
        lmp.command("create_atoms 0 random 1 %s NULL mol methane %s" % (randnum(), randnum()))
        print "test"
        lmp.group("del type 3")
        print "test"
        run(0)
        if rank == 0:
            bolz.append(math.exp(-(lmp.eval("pe")-E_step)/(natoms*kb*T)))
        lmp.delete_atoms("group del")
    # Calculates step - based quantities.
    if rank == 0:
        Bi_step = np.sum(bolz) * vol / ninsert
        MUex_step = np.sum(bolz) / ninsert
        av_vol += vol
        av_rho += rho
        Bi += Bi_step
        MUex += MUex_step
        # Print out to a log file.
        log.write("%s %s" % ((i+1),-kb*T*math.log(MUex/(i+1))))
if rank == 0:    
    # Calcualate the average values 
    Bi /= nsteps
    MUex /= nsteps
    av_vol /= nsteps
    av_rho /= nsteps

    # Convert to 1/(kg*bar)
    nistconv = 0.01438329989 
    MUinf = -kb*T/av_vol*math.log(Bi)
    MUex = -kb*T*math.log(MUex)
    H = av_rho*kb*T*math.exp(MUinf/(kb*T))*nistconv
    print MUinf,MUex,H 
