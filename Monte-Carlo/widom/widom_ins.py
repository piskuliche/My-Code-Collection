import numpy as np
import sys, random, math
from mpi4py import MPI

"""
This is a code for doing widom insertions inside LAMMPS using the python wrapper lammps.py

Usage: python widom_ins.py infile ninsert nsteps molfile

"""

def run(steps):
    lmp.command("run %s" % int(steps))

def grabpe():
    pe = lmp.extract_compute("thermo_pe",0,0)/natoms
    return pe

def grabvol():
    vol = (lmp.extract_global("boxxhi",1)-lmp.extract_global("boxxlo",1))**3
    return vol


# Define Units
kb   = 0.0019872041 # kcal/(mol*K)
T    = 314.1 # K
beta = (1.0/(kb*T))

"""
convfactor takes:
    1       kcal
    -   *   ---- => kbar
    A^3      mol
"""
convfactor = 69.47511748

# MPI Section
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Argument Parser
if len(sys.argv) != 5:
    print("Usage: python widom_ins.py infile ninsert nsteps molfile")
    sys.exit(1)
infile  = sys.argv[1]
ninsert = int(sys.argv[2])
nsteps  = int(sys.argv[3])
molfile = sys.argv[4]

# Generate Random Number Seed
random.seed(27413)

# Define Random Number Generator
def randnum():
    rand = ''.join(random.sample("0123456789",6))
    return rand

# Import LAMMPS and initialize it
from lammps import lammps
lmp = lammps()

# Read LAMMPS inputfile line by line
lines = open(infile, 'r').readlines()
for line in lines: lmp.command(line)

run(0)
lmp.command("variable e equal pe")

# Read initial information
natoms = lmp.extract_global("natoms", 0)
inite = grabpe()
vol = grabvol()

# Print Basic info to screen.
if rank == 0:
    print("%s atoms in the system" % natoms)
    print("%s is the initial energy (kcal/mol)" % inite)

lmp.command("molecule ins %s toff 2" % molfile)
run(100000)

# Initialize
vol = []
nrho = []
Bi = []
# Begin Step Loop
for i in range(nsteps):
    # Zeros the energy vlaues
    bolz = []
    # Move to new config
    run(1000)
    E_h2o = grabpe()
    vol.append(grabvol())
    nrho.append(natoms/vol[i])
    # Loops over insertion attempts.
    for j in range(ninsert):
        lmp.command("create_atoms 0 random 1 %s NULL mol ins %s" % (randnum(),randnum()))
        lmp.command("group del type 3")
        run(0)
        ins_e = grabpe()
        del_e = ins_e - E_h2o
        bolz.append(math.exp(-beta*del_e))
        lmp.command("delete_atoms group del")
    # Calculates Quantities for the step.
    Bi.append(np.average(bolz)*vol[i])

# Open the Logfile
logfile = "log.widom"
log = open(logfile,'w')
if rank == 0:
    log.write("Step MuEx H\n")

av_vol = np.average(vol)
av_nrho = np.average(nrho)
for i in range(nsteps):
    MuEx = -kb*T*math.log(Bi[i]/av_vol)
    H    = av_nrho/beta*math.exp(beta*MuEx)
    if rank == 0:
        log.write("%s %s %s\n" % ((i+1),MuEx, H*convfactor))
        log.flush()
 
MuEx = -kb*T*math.log(np.average(Bi)/av_vol)
H = av_nrho/beta*math.exp(beta*MuEx)
log.close()
if rank == 0:
    print("MuEx(kcal/mol) %s" % (MuEx))
    print("Henry(kbar) %s" % (H*convfactor))


    
