#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os, json
import molecular_species as molspec

"""
This is a python program which can build a box of generally defined molecules.
These molecules cane be included in molecular_species.py
"""

if not os.path.exists('./tmp'):
        os.makedirs('./tmp')

# Read the NML file
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys, and typcheck values
defaults        = {"nspec":[0], "tspec":["bmim"],"blength":40.0,"num_components":1}
    
# Set parameters to input values or defaults
num_spec       = nml["nspec"]          if "nspec"          in nml else defaults["nspec"]
typ_spec       = nml["tspec"]          if "tspec"          in nml else defaults["tspec"]
blength         = nml["blength"]        if "blength"        in nml else defaults["blength"]
num_components  = nml["num_components"] if "num_components" in nml else defaults["num_components"]

# Write out parameters
print("There are %s total components" % num_components)
for i in range(num_components):
    print("Species %s has %s %s molecules." % (i, num_spec[i],typ_spec[i]))
print("box length = %s" % blength)

# Define functions
def gen_packmol():
    L = blength-2.0
    pmolfile = 'system.pmol'
    pm = open(pmolfile, 'w')
    pm.write("tolerance 2.0\n")
    pm.write("filetype xyz\n")
    pm.write("output system.xyz\n")
    pm.write("\n")
    pm.close()
    NCHAR =[]
    NTYPS =[]
    ATMS = []
    BNDS = []
    ANGS = []
    DIHS = []
    IMPS = []
    for i in range(num_components):
        nchar,ntyps,atms,bnds,angs,dihs,imps=molspec.molecule(typ_spec[i], num_spec[i], blength)
        NCHAR.append(nchar)
        NTYPS.append(ntyps)
        ATMS.append(atms)
        BNDS.append(bnds)
        ANGS.append(angs)
        DIHS.append(dihs)
        IMPS.append(imps)
    return NCHAR, NTYPS, ATMS, BNDS, ANGS, DIHS, IMPS

def write_header():
    datafile = "data.lmps"
    dat      = open(datafile, 'w')

    dat.write("Data file with components " )
    for i in range(num_components):
        dat.write("%s " % typ_spec[i])
    dat.write("\n\n")
    nATMS, nBNDS, nANGS, nDIHS, nIMPS = 0,0,0,0,0
    natyps, nbntyps, nantyps, ndityps, nimtyps = 0,0,0,0,0
    for i in range(num_components):
        nATMS += NCHAR[i][0]*num_spec[i] 
        natyps += NTYPS[i][0]
        nBNDS += NCHAR[i][1]*num_spec[i]
        nbntyps += NTYPS[i][1]
        nANGS += NCHAR[i][2]*num_spec[i]
        nantyps += NTYPS[i][2] 
        nDIHS += NCHAR[i][3]*num_spec[i]
        ndityps += NTYPS[i][3] 
        nIMPS += NCHAR[i][4]*num_spec[i]
        nimtyps += NTYPS[i][4] 
    dat.write("%s atoms\n" % nATMS)
    dat.write("%s bonds\n" % nBNDS)
    dat.write("%s angles\n" % nANGS)
    dat.write("%s dihedrals\n" % nDIHS)
    dat.write("%s impropers\n" % nIMPS)
    dat.write("\n")
    dat.write("%s atom types\n" % np.sum(natyps))
    dat.write("%s bond types\n" % np.sum(nbntyps))
    dat.write("%s angle types\n" % np.sum(nantyps))
    dat.write("%s dihedral types\n" % np.sum(ndityps))
    dat.write("%s improper types\n" % np.sum(nimtyps))
    dat.write("\n")
    dat.write("0.0   %s xlo xhi\n" % blength)
    dat.write("0.0   %s ylo yhi\n" % blength)
    dat.write("0.0   %s zlo zhi\n" % blength)
    dat.write("\n")
    # Write the masses section of the data file
    dat.write("Masses\n")
    dat.write("\n")
    count = 1
    # Species 1 Masses
    for i in range(num_components):
        a, indices = np.unique(ATMS[i]["atype"],return_index=True)
        for j in indices:
            dat.write("%s %s\n" % (count, ATMS[i]["mass"][j]))
            count += 1
    dat.write("\n")
    dat.close()
    return

def write_atoms():
    datafile = "data.lmps"
    dat      = open(datafile, 'a')
    dat.write("Atoms\n")
    dat.write("\n")
    shiftatoms = 0
    shiftmax = 0
    molcount = 0
    for i in range(num_components):
        for mol in range(num_spec[i]):
            for atom in range(NCHAR[i][0]):
                mindex = molcount + 1
                aindex = mol*NCHAR[i][0] + 1 + atom + shiftatoms
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[i]["atype"][atom]+shiftmax, ATMS[i]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            molcount += 1
        shiftatoms += num_spec[i]*NCHAR[i][0]
        shiftmax += max(ATMS[i]["atype"])

    dat.write("\n")
    dat.close()
    return

def write_bonds():
    datafile = "data.lmps"
    dat      = open(datafile, 'a')
    dat.write("Bonds\n")
    dat.write("\n")
    count, typstart, shiftbnds, setbndsht = 1, 0, 0, 0
    # Writes the bonding section
    for i in range(num_components):
        for mol in range(num_spec[i]):
            for bnd in BNDS[i]["name"]:
                shiftbnds = NCHAR[i][0]*mol + 1 + setbndsht
                dat.write("%s %s %s %s\n" % (count,typstart+BNDS[i][str(bnd)][0],BNDS[i][str(bnd)][3]+shiftbnds,BNDS[i][str(bnd)][4]+shiftbnds))
                count += 1
        typstart += NTYPS[i][1]
        setbndsht += NCHAR[i][0]*num_spec[i]
    dat.write("\n")
    dat.close()

def write_angles():
    datafile = "data.lmps"
    dat      = open(datafile, 'a')
    dat.write("Angles\n")
    dat.write("\n")
    count, typstart, shiftangs, setangsht = 1, 0, 0, 0
    for i in range(num_components):
        for mol in range(num_spec[i]):
            for ang in ANGS[i]['name']:
                shiftangs = NCHAR[i][0]*mol + 1 + setangsht
                dat.write("%s %s %s %s %s\n" % (count, typstart+ANGS[i][str(ang)][0], ANGS[i][str(ang)][3]+shiftangs,ANGS[i][str(ang)][4]+shiftangs,ANGS[i][str(ang)][5]+shiftangs))
                count += 1
        typstart += NTYPS[i][2]
        setangsht += NCHAR[i][0]*num_spec[i]
    dat.write("\n")
    dat.close()

def write_dihedrals():
    datafile = "data.lmps"
    dat      = open(datafile, 'a')
    dat.write("Dihedrals\n")
    dat.write("\n")
    count, typstart, shiftdihs, setdihsht = 1, 0, 0, 0
    for i in range(num_components):
        for mol in range(num_spec[i]):
            for dih in DIHS[i]['name']:
                shiftdihs = NCHAR[i][0]*mol + 1 + setdihsht
                dat.write("%s %s %s %s %s %s\n" % (count, typstart+DIHS[i][str(dih)][0], DIHS[i][str(dih)][4]+shiftdihs,DIHS[i][str(dih)][5]+shiftdihs,DIHS[i][str(dih)][6]+shiftdihs,DIHS[i][str(dih)][7]+shiftdihs))
                count += 1
        typstart += NTYPS[i][3]
        setdihsht += NCHAR[i][0]*num_spec[i]
    dat.write("\n")
    dat.close()

def write_impropers():
    datafile = "data.lmps"
    dat      = open(datafile, 'a')
    dat.write("Impropers\n")
    dat.write("\n")
    count, typstart, shiftimps, setimpsht = 1, 0, 0, 0
    for i in range(num_components):
        for mol in range(num_spec[i]):
            for imp in IMPS[i]['name']:
                shiftimps = NCHAR[i][0]*mol + 1 + setimpsht
                dat.write("%s %s %s %s %s %s\n" % (count, typstart+IMPS[i][str(imp)][0], IMPS[i][str(imp)][3]+shiftimps,IMPS[i][str(imp)][4]+shiftimps,IMPS[i][str(imp)][5]+shiftimps,IMPS[i][str(imp)][6]+shiftimps))
                count += 1
        typstart += NTYPS[i][4]
        setimpsht += NCHAR[i][0]*num_spec[i]
    dat.write("\n") 
    dat.write("\n")
    dat.close()

def write_lammps_paircoeffs():
    lmpsfile = 'lmps.paircoeffs'
    lmps = open(lmpsfile, 'w')
    if 'tip4p' not in typ_spec:
        lmps.write("pair_style lj/cut/coul/long 12.5 12.5\n")
    else:
        lmps.write("#!!!! NEED TO USE A HYBRID PAIRSTYLE -> need to edit pair_coeff commands accordingly\n")
        lmps.write("pair_style hybrid lj/cut/tip4p/long 12.5 12.5 lj/cut/coul/long 12.5 12.5\n")
    start = 0
    for i in range(num_components):
        lmps.write("# Pair Coeffs Species %s\n" % i)
        a, indices = np.unique(ATMS[i]["atype"], return_index = True)
        for j in indices:
            lmps.write("pair_coeff %s %s %9.3f %9.3f\n" % (ATMS[i]["atype"][j]+start,ATMS[i]["atype"][j]+start, ATMS[i]["kjeps"][j]/4.184, ATMS[i]["rmin"][j]/(2**(1/6.))))
        start += max(ATMS[i]["atype"])
    lmps.close()
    return

def write_lammps_bondcoeffs():
    lmpsfile = 'lmps.bondcoeffs'
    lmps = open(lmpsfile, 'w')
    lmps.write("bond_style harmonic\n")
    start = 0
    usedtyps=[]
    for i in range(num_components):
        lmps.write("# Bond Coeffs Species %s\n" % i)
        for bond in BNDS[i]['name']:
            if BNDS[i][bond][0]+start not in usedtyps:
                lmps.write("bond_coeff %s %9.3f %9.3f\n" % (BNDS[i][bond][0]+start, BNDS[i][bond][1], BNDS[i][bond][2]))
                usedtyps.append(BNDS[i][bond][0]+start)
        start = max(usedtyps)
    lmps.close()
    return

def write_lammps_anglecoeffs():
    lmpsfile = 'lmps.anglecoeffs'
    lmps = open(lmpsfile, 'w')
    lmps.write("angle_style harmonic\n")
    start = 0
    usedtyps=[]
    for i in range(num_components):
        lmps.write("# Angle Coeffs Species %s\n" % i)
        for angle in ANGS[i]['name']:
            if ANGS[i][angle][0]+start not in usedtyps:
                lmps.write("angle_coeff %s %9.3f %9.3f\n" % (ANGS[i][angle][0]+start, ANGS[i][angle][1], ANGS[i][angle][2]))
                usedtyps.append(ANGS[i][angle][0]+start)
        start = max(usedtyps)
    lmps.close()
    return

def write_lammps_dihedralcoeffs():
    lmpsfile = 'lmps.dihedralcoeffs'
    lmps = open(lmpsfile, 'w')
    lmps.write("dihedral_style charmm\n")
    start = 0
    usedtyps=[]
    for i in range(num_components):
        lmps.write("# Dihedral Coeffs Species %s\n" % i)
        for dihedral in DIHS[i]['name']:
            if DIHS[i][dihedral][0]+start not in usedtyps:
                lmps.write("dihedral_coeff %s %9.3f %s %9.3f\n" % (DIHS[i][dihedral][0]+start, DIHS[i][dihedral][1], DIHS[i][dihedral][2], DIHS[i][dihedral][3]))
                usedtyps.append(DIHS[i][dihedral][0]+start)
        start = max(usedtyps)
    lmps.close()
    return

def write_lammps_impropercoeffs():
    lmpsfile = 'lmps.impropercoeffs'
    lmps = open(lmpsfile, 'w')
    lmps.write("improper_style harmonic\n")
    start = 0
    usedtyps=[]
    for i in range(num_components):
        lmps.write("# Improper Coeffs Species %s\n" % i)
        for improper in IMPS[i]['name']:
            if IMPS[i][improper][0]+start not in usedtyps:
                lmps.write("improper_coeff %s %9.3f %9.3f\n" % (IMPS[i][improper][0]+start, IMPS[i][improper][1], IMPS[i][improper][2]))
                usedtyps.append(IMPS[i][improper][0]+start)
        start = max(usedtyps)
    lmps.close()
    return







NCHAR, NTYPS, ATMS, BNDS, ANGS, DIHS, IMPS = gen_packmol()
#subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < system.pmol"], shell=True)
x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
write_header()
write_atoms()
write_bonds()
write_angles()
write_dihedrals()
write_impropers()
write_lammps_paircoeffs()
write_lammps_bondcoeffs()
write_lammps_anglecoeffs()
write_lammps_dihedralcoeffs()
write_lammps_impropercoeffs()
