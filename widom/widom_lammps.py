import numpy as np
import sys,random,math
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print('My rank is ',rank)

kb=0.0019872041 #kcal/molK
T=298.15 # K
ninsert = 200
nsteps = 100
random.seed(27413)

def randnum():
    rand = val = ''.join(random.sample("0123456789", 6))
    return rand

# Read Command Line Arguments
argv = sys.argv
if len(argv) != 2:
    print("Syntax: mc.py in.mc")
    sys.exit()

infile=sys.argv[1]

from lammps import lammps
lmp = lammps()

# read in input file

lines = open(infile,'r').readlines()
for line in lines: lmp.command(line)

lmp.command("run 0")

lmp.command("variable e equal pe")

# calculate the energy of the system
lmp.command("run 0")

natoms = lmp.extract_global("natoms", 0)
inite = lmp.extract_compute("thermo_pe", 0,0,) / natoms
vol = lmp.extract_global("boxxhi",1)**3


if rank==0:
    print("%s atoms in the system" % natoms)
    print("%s is the initial energy (kcal/mol)" % inite)

lmp.command("molecule methane ch4.txt toff 2")
lmp.command("run 100000")

# Array Parameters
"""
         bolz : e^-du/kbT array over ninsert
      Bi_step : step value of V*e^-du/kbT
           Bi : sum over steps of V*e^-du/kbT
    MUex_step : step value of e^-du/kbT
         MUex : Excess chemical potential 
       E_step : The energy of the step (before insertions)
"""

Bi_step, Bi, MUex_step, MUex, MUinf = 0.0, 0.0, 0.0, 0.0, 0.0
av_vol, av_rho = 0.0, 0.0
E_step = 0.0

# Open logfile
logfile = "log.widom"
log = open(logfile, 'w')
log.write("Step MUinf MUex\n")
# Begin step loop
for i in range(nsteps):
    # Zeros the energy values.
    bolz = []
    Bi_step, MUex_step = 0.0, 0.0
    # Runs 1000 steps, computes pe, grabs volume, grabs density
    lmp.command("run 1000")
    E_step=lmp.extract_compute("thermo_pe",0,0)
    vol = lmp.extract_global("boxxhi",1)**3
    rho = lmp.get_thermo("density")
    # Loop over insertion attempts.
    for i in range(ninsert):
        lmp.command("create_atoms 0 random 1 %s NULL mol methane %s" % (randnum(),randnum()))
        lmp.command("group del type 3")
        lmp.command("run 0")
        bolz.append(math.exp(-(lmp.extract_compute("thermo_pe",0,0)-E_step)/(natoms*kb*T)))
        lmp.command("delete_atoms group del")
    # Calculates step - based quantities.
    Bi_step = np.sum(bolz) * vol / ninsert
    MUex_step = np.sum(bolz) / ninsert
    av_vol += vol
    av_rho += rho
    Bi += Bi_step
    MUex += MUex_step
    # Print out to a log file.
    log.write("%s %s" % ((i+1),-kb*T*math.log(MUex/(i+1))))
    
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
