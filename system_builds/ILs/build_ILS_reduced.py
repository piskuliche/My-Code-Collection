#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os, json
import molecular_species as molspec

"""
This is a python program which can build a box of ionic liquids for use in lammps simulations.
"""

# Read the NML file
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys, and typcheck values
defaults        = {"nspec1":0, "nspec2":0, "nspec3":0, "tspec1":"bmim", "tspec2":"pf6", "tspec3":"R32","blength":40.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]),key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))
    
# Set parameters to input values or defaults
num_spec1       = nml["nspec1"]         if "nspec1"         in nml else defaults["nspec1"]
num_spec2      = nml["nspec2"]        if "nspec2"        in nml else defaults["nspec2"]
num_spec3      = nml["nspec3"]        if "nspec3"        in nml else defaults["nspec3"]
typ_spec1       = nml["tspec1"]      if "tspec1"      in nml else defaults["tspec1"]
typ_spec2      = nml["tspec2"]     if "tspec2"     in nml else defaults["tspec2"]
typ_spec3      = nml["tspec3"]     if "tspec3"     in nml else defaults["tspec3"]
blength         = nml["blength"]     if "blength"     in nml else defaults["blength"]

# Write out parameters
print("num_spec1 = %s of type %s" % (num_spec1,typ_spec1))
print("num_spec2 = %s of type %s" % (num_spec2,typ_spec2))
print("num_spec3 = %s of type %s" % (num_spec3,typ_spec3))
print("box length = %s" % blength)

# Define functions
def gen_packmol():
    L = blength-2.0
    pmolfile = 'IL_system.pmol'
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
    if num_spec1 != 0:
        nchar,ntyps,atms,bnds,angs,dihs,imps=molspec.molecule(typ_spec1, num_spec1, blength)
        NCHAR.append(nchar)
        NTYPS.append(ntyps)
        ATMS.append(atms)
        BNDS.append(bnds)
        ANGS.append(angs)
        DIHS.append(dihs)
        IMPS.append(imps)
    if num_spec2 != 0:
        nchar,ntyps,atms,bnds,angs,dihs,imps=molspec.molecule(typ_spec2, num_spec2, blength)
        NCHAR.append(nchar)
        NTYPS.append(ntyps)
        ATMS.append(atms)
        BNDS.append(bnds)
        ANGS.append(angs)
        DIHS.append(dihs)
        IMPS.append(imps)
    if num_spec3 != 0:
        nchar,ntyps,atms,bnds,angs,dihs,imps=molspec.molecule(typ_spec3, num_spec2, blength)
        NCHAR.append(nchar)
        NTYPS.append(ntyps)
        ATMS.append(atms)
        BNDS.append(bnds)
        ANGS.append(angs)
        DIHS.append(dihs)
        IMPS.append(imps)
    return NCHAR, NTYPS, ATMS, BNDS, ANGS, DIHS, IMPS

def write_header():
    datafile = "data.IL"
    dat      = open(datafile, 'w')
    dat.write("Ionic Liquid Data File for %s %s and %s\n" % (typ_spec1, typ_spec2, typ_spec3))
    dat.write("\n")
    nATMS = NCHAR[0][0]*num_spec2 + NCHAR[1][0]*num_spec1 + NCHAR[2][0]*num_spec3
    natyps = NTYPS[0][0] + NTYPS[1][0] + NTYPS[2][0]
    nBNDS = NCHAR[0][1]*num_spec2 + NCHAR[1][1]*num_spec1 + NCHAR[2][1]*num_spec3
    nbtyps = NTYPS[0][1] + NTYPS[1][1] + NTYPS[2][1]
    nANGS = NCHAR[0][2]*num_spec2 + NCHAR[1][2]*num_spec1 + NCHAR[2][2]*num_spec3
    natyps = NTYPS[0][2] + NTYPS[1][2] + NTYPS[2][2]
    nDIHS = NCHAR[0][3]*num_spec2 + NCHAR[1][3]*num_spec1 + NCHAR[2][3]*num_spec3
    ndtyps = NTYPS[0][3] + NTYPS[1][3] + NTYPS[2][3]
    nIMPS = NCHAR[0][4]*num_spec2 + NCHAR[1][4]*num_spec1 + NCHAR[2][4]*num_spec3
    nityps = NTYPS[0][4] + NTYPS[1][4] + NTYPS[2][4]
    dat.write("%s atoms\n" % nATMS)
    dat.write("%s bonds\n" % nBNDS)
    dat.write("%s angles\n" % nANGS)
    dat.write("%s dihedrals\n" % nDIHS)
    dat.write("%s impropers\n" % nIMPS)
    dat.write("\n")
    dat.write("%s atom types\n" % np.sum(natyps))
    dat.write("%s bond types\n" % np.sum(nbtyps))
    dat.write("%s angle types\n" % np.sum(natyps))
    dat.write("%s dihedral types\n" % np.sum(ndtyps))
    dat.write("%s improper types\n" % np.sum(nityps))
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
    a, indices = np.unique(ATMS[0]["atype"],return_index=True)
    for j in indices:
        print(j)
        dat.write("%s %s\n" % (count, ATMS[0]["mass"][j]))
        count += 1
    a, indices = np.unique(ATMS[1]["atype"],return_index=True)
    for j in indices:
        print(j)
        dat.write("%s %s\n" % (count, ATMS[1]["mass"][j]))
        count += 1
    a, indices = np.unique(ATMS[2]["atype"],return_index=True)
    for j in indices:
        print(j)
        dat.write("%s %s\n" % (count, ATMS[2]["mass"][j]))
        count += 1
    dat.close()
    return

def write_atoms():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Atoms\n")
    dat.write("\n")
    for mol in range(num_spec1+num_spec2+num_spec3):
        if mol < num_spec1:
            moltype = 0
        elif mol >= num_spec1 and mol < num_spec1+num_spec2:
            moltype = 1
        else:
            moltype = 2
        for atom in range(NCHAR[moltype][0]):
            mindex = mol + 1
            if moltype == 0:
                aindex   = mol*NCHAR[moltype][0] + 1 + atom
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[[moltype]["atype"][atom], ATMS[moltype]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            elif moltype == 1:
                aindex = (mol-num_spec1)*NCHAR[moltype][0] + 1 + atom + num_spec1*NCHAR[moltype-1][0]
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[moltype]["atype"][atom]+max(ATMS[moltype-1]["atype"][0]), ATMS[moltype]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            else:
                aindex = (mol-num_spec2-num_spec1)*NCHAR[moltype][0] + 1 + atom + num_spec1*NCHAR[moltype-2][0] + num_spec2*NCHAR[moltype-1][1]
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[moltype]["atype"][atom]+max(ATMS[moltype-1]["atype"][0])+max(ATMS[moltype-2]["atype"][0]), ATMS[moltype]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
    dat.write("\n")
    dat.close()
    return

def write_bonds():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Bonds\n")
    dat.write("\n")
   for mol in range(num_spec2):
        aindex   = mol*natms[0] + 1 
        bndstart = mol*nbnds[0]
        typstart = 0
        line, line2, line3, line4 = spec2_bonding(aindex-1, bndstart,0,0,0, typstart,0,0,0)
        for i in range(nbnds[0]):
            dat.write(line[i])
    for mol in range(num_spec2, num_spec2+num_spec1):
        aindex   = (mol-num_spec2)*natms[1] + 1 + num_spec2*natms[0]
        bndstart = (mol-num_spec2)*nbnds[1] + num_spec2*nbnds[0]
        typstart = nbtyps[0]
        line, line2= spec1_bonding(aindex-1, bndstart,0, typstart,0)
        for i in range(nbnds[1]):
            dat.write(line[i])
    for mol in range(num_spec2+num_spec1, num_spec2+num_spec1+num_spec3):
        aindex   = (mol-num_spec2-num_spec1)*natms[2] + 1 + num_spec2*natms[0] + num_spec1*natms[1]
        bndstart = (mol-num_spec2-num_spec1)*nbnds[2] + num_spec2*nbnds[0] + num_spec1*nbnds[1]
        typstart = nbtyps[0] + nbtyps[1]
        line, line2 = spec3_bonding(aindex-1, bndstart,0, typstart,0)
        for i in range(nbnds[2]):
            dat.write(line[i])
    dat.write("\n")
    dat.close()

def write_angles():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Angles\n")
    dat.write("\n")
    for mol in range(num_spec2):
        aindex = mol*natms[0] + 1
        bndstart, angstart = mol*nbnds[0], mol*nangs[0]
        typstart = 0
        angtypstart = 0
        line, line2,line3,line4 = spec2_bonding(aindex-1, bndstart, angstart,0,0, typstart, angtypstart,0,0)
        for i in range(nangs[0]):
            dat.write(line2[i])
    for mol in range(num_spec2, num_spec2+num_spec1):
        aindex   = (mol-num_spec2)*natms[1] + 1 + num_spec2*natms[0]
        bndstart, angstart = (mol-num_spec2)*nbnds[1] + num_spec2*nbnds[0], (mol-num_spec2)*nangs[1] + num_spec2*nangs[0]
        typstart = nbtyps[0]
        angtypstart = nangtyps[0]
        print(angtypstart)
        line,line2 = spec1_bonding(aindex-1, bndstart, angstart, typstart, angtypstart)
        for i in range(nangs[1]):
            dat.write(line2[i])
    for mol in range(num_spec2+num_spec1, num_spec2+num_spec1+num_spec3):
        aindex   = (mol-num_spec2-num_spec1)*natms[2] + 1 + num_spec2*natms[0] + num_spec1*natms[1]
        bndstart,angstart = (mol-num_spec2-num_spec1)*nbnds[2] + num_spec2*nbnds[0] + num_spec1*nbnds[1],(mol-num_spec2-num_spec1)*nangs[2] + num_spec2*nangs[0] + num_spec1*nangs[1]
        typstart = nbtyps[0] + nbtyps[1]
        angtypstart = nangtyps[0]+nangtyps[1]
        line, line2 = spec3_bonding(aindex-1, bndstart,angstart, typstart,angtypstart)
        for i in range(nangs[2]):
            dat.write(line2[i])
    dat.write("\n")
    dat.close()

def write_dihedrals():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Dihedrals\n")
    dat.write("\n")
    for mol in range(num_spec2):
        aindex = mol*natms[0] + 1
        bndstart, angstart, dihstart= mol*nbnds[0], mol*nangs[0], mol*ndihs[0]
        typstart = 0
        angtypstart = 0
        dihtypstart = 0
        line, line2, line3, line4 = spec2_bonding(aindex-1, bndstart, angstart,dihstart,0, typstart, angtypstart,dihtypstart,0)
        for i in range(ndihs[0]):
            dat.write(line3[i])
    dat.write("\n")
    dat.close()

def write_impropers():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Impropers\n")
    dat.write("\n")
    for mol in range(num_spec2):
        aindex = mol*natms[0] + 1
        bndstart, angstart, dihstart, impstart = mol*nbnds[0], mol*nangs[0], mol*ndihs[0], mol*nimps[0]
        typstart = 0
        angtypstart = 0
        dihtypstart = 0
        imptypstart = 0
        line, line2, line3, line4 = spec2_bonding(aindex-1, bndstart, angstart,dihstart,impstart, typstart, angtypstart,dihtypstart,imptypstart)
        for i in range(nimps[0]):
            dat.write(line4[i])
    dat.write("\n")
    dat.close()




NCHAR, NTYPS, ATMS, BNDS, ANGS, DIHS, IMPS = gen_packmol()
print(NCHAR)
print(ATMS[0])
print(ATMS[1])
print(ATMS[2])
#subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < IL_system.pmol"], shell=True)
#x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
write_header()
#write_atoms()
#write_bonds()
#write_angles()
#write_dihedrals()
#write_impropers()
