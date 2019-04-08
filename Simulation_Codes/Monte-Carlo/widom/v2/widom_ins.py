import numpy as np
import sys, random, math
from mpi4py import MPI

"""
This is a code for doing widom insertions inside LAMMPS using the python wrapper lammps.py

Usage: python widom_ins.py infile ninsert nsteps molfile

"""

#-------------------------------------------------------------------------------
"""
This includes all the fucntion definitions.
"""
#-------------------------------------------------------------------------------

# Define Random Number Generator
def randnum():
    rand = ''.join(random.sample("0123456789",6))
    return rand

# Makes Lammps run steps number of timesteps.
def run(steps):
    lmp.command("run %s" % int(steps))

# Grabs the energy from lammps
def grabpe():
    pe = lmp.extract_compute("thermo_pe",0,0)/natoms
    return pe

# Grabs the volume from lammps
def grabvol():
    vol = (lmp.extract_global("boxxhi",1)-lmp.extract_global("boxxlo",1))**3
    return vol

# Writes an example input file
def write_pyf():
    pyf = open('input.template','w')
    pyf.write("# Number of Configurations\n\n")
    pyf.write("# Number of Insertions per Step\n\n")
    pyf.write("# Molecule File to Insert\n\n")
    pyf.write("# Temperature (K)\n\n")
    pyf.write("# Pressure (bar)\n\n")
    pyf.close() 

# Reads the input file
def read_pyf(file):
    inp_params = np.genfromtxt(file, comments="#", usecols=0)
    print("Number of Configurations: %s" % inp_params[0])
    print("Number of Insertions:     %s" % inp_params[1])
    print("Molecule File to Insert:  %s" % inp_params[2])
    print("Temperature (K): %s" % inp_params[3])
    print("Pressure (bar)" % inp_params[4])
    return inp_params


#-------------------------------------------------------------------------------
"""
This is the header information for the program. 

It also reads in all of the file input parameters.
"""
#-------------------------------------------------------------------------------
# Define Units
kb   = 0.0019872041 # kcal/(mol*K)
T    = 298.15 # K
P    = 1 # bar
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
if len(sys.argv) != 3:
    print("Usage: python widom_ins.py infile pyfile")
    print("If you need an input file exampe, set pyfile to FALSE")
    sys.exit(1)
infile  = sys.argv[1]
pyfile  = sys.argv[2]

if pyfile == "FALSE":
    print("Pyfile Generated")
    write_pyf()
    sys.exit(0)

# Read input parameters
pyin = read_pyf(pyfile)
nsteps  = pyin[0]
nins    = pyin[1]
molfile = pyin[2]
T       = pyin[3]
P       = pyin[4]


# Generate Random Number Seed
random.seed(27413)

#-------------------------------------------------------------------------------
"""
This section initializes lammps and runs the commands included in the initial script.

It runs 100 picoseconds of equilibration in whichever ensemble is originally initialized.
"""
#-------------------------------------------------------------------------------
# Initialize Arrays
vol = []
nrho = []
Bi = []

# Import LAMMPS and initialize it
from lammps import lammps
lmp = lammps()

# Read LAMMPS inputfile line by line
lines = open(infile, 'r').readlines()
for line in lines: lmp.command(line)

# Calculates initial energies
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
#-------------------------------------------------------------------------------
"""
This section handles the actual loop over configurations and insertions. 
"""
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
"""
SECTION: OUTPUT

This section does all the outputting for the simulation.
Outputs: log.widom
    Stores widom information: step MuEx H Mu
"""
#-------------------------------------------------------------------------------
# Open the Logfile
logfile = "log.widom"
log = open(logfile,'w')
if rank == 0:
    log.write("Step MuEx H Mu\n")

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


    
