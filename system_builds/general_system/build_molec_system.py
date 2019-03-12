#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os, json
import molecular_species as molspec

"""
This is a python program which can build a box of generally defined molecules.
These molecules cane be included in molecular_species.py
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
        if val != nml["tspec1"] and val != nml["tspec2"] and val != nml["tspec3"]:
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
    nbntyps = NTYPS[0][1] + NTYPS[1][1] + NTYPS[2][1]
    nANGS = NCHAR[0][2]*num_spec2 + NCHAR[1][2]*num_spec1 + NCHAR[2][2]*num_spec3
    nantyps = NTYPS[0][2] + NTYPS[1][2] + NTYPS[2][2]
    nDIHS = NCHAR[0][3]*num_spec2 + NCHAR[1][3]*num_spec1 + NCHAR[2][3]*num_spec3
    ndityps = NTYPS[0][3] + NTYPS[1][3] + NTYPS[2][3]
    nIMPS = NCHAR[0][4]*num_spec2 + NCHAR[1][4]*num_spec1 + NCHAR[2][4]*num_spec3
    nimtyps = NTYPS[0][4] + NTYPS[1][4] + NTYPS[2][4]
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
    a, indices = np.unique(ATMS[0]["atype"],return_index=True)
    for j in indices:
        dat.write("%s %s\n" % (count, ATMS[0]["mass"][j]))
        count += 1
    a, indices = np.unique(ATMS[1]["atype"],return_index=True)
    for j in indices:
        dat.write("%s %s\n" % (count, ATMS[1]["mass"][j]))
        count += 1
    a, indices = np.unique(ATMS[2]["atype"],return_index=True)
    for j in indices:
        dat.write("%s %s\n" % (count, ATMS[2]["mass"][j]))
        count += 1
    dat.write("\n")
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
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[moltype]["atype"][atom], ATMS[moltype]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            elif moltype == 1:
                aindex = (mol-num_spec1)*NCHAR[moltype][0] + 1 + atom + num_spec1*NCHAR[moltype-1][0]
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[moltype]["atype"][atom]+max(ATMS[moltype-1]["atype"]), ATMS[moltype]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            else:
                aindex = (mol-num_spec2-num_spec1)*NCHAR[moltype][0] + 1 + atom + num_spec1*NCHAR[moltype-2][0] + num_spec2*NCHAR[moltype-1][0]
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, ATMS[moltype]["atype"][atom]+max(ATMS[moltype-1]["atype"])+max(ATMS[moltype-2]["atype"]), ATMS[moltype]["q"][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
    dat.write("\n")
    dat.close()
    return

def write_bonds():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Bonds\n")
    dat.write("\n")
    count = 1
    for mol in range(num_spec1):
        typstart = 0
        for i in BNDS[0]['name']:
            shift = NCHAR[0][0]*mol+1
            dat.write("%s %s %s %s\n" % (count,typstart+BNDS[0][str(i)][0],BNDS[0][str(i)][3]+shift,BNDS[0][str(i)][4]+shift))
            count += 1
    for mol in range(num_spec2, num_spec2+num_spec1):
        typstart = NTYPS[0][2]
        for i in BNDS[1]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*(mol-num_spec1) + 1
            dat.write("%s %s %s %s\n" % (count,typstart+BNDS[1][str(i)][0],BNDS[1][str(i)][3]+shift,BNDS[1][str(i)][4]+shift))
            count += 1
    for mol in range(num_spec1+num_spec2, num_spec1+num_spec2+num_spec3):
        typstart = NTYPS[0][2] + NTYPS[1][2]
        for i in BNDS[2]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*num_spec2 + NCHAR[2][0]*(mol-num_spec1-num_spec2) +1
            dat.write("%s %s %s %s\n" % (count,typstart+BNDS[2][str(i)][0],BNDS[2][str(i)][3]+shift,BNDS[2][str(i)][4]+shift))
            count += 1
    dat.write("\n")
    dat.close()

def write_angles():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Angles\n")
    dat.write("\n")
    count = 1
    for mol in range(num_spec1):
        typstart = 0
        for i in ANGS[0]['name']:
            shift = NCHAR[0][0]*mol + 1
            dat.write("%s %s %s %s %s\n" % (count, typstart+ANGS[0][str(i)][0], ANGS[0][str(i)][3]+shift,ANGS[0][str(i)][4]+shift,ANGS[0][str(i)][5]+shift))
            count += 1
    for mol in range(num_spec2, num_spec2+num_spec1):
        typstart = NTYPS[0][2]
        for i in ANGS[1]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*(mol-num_spec1) + 1
            dat.write("%s %s %s %s %s\n" % (count, typstart+ANGS[1][str(i)][0], ANGS[1][str(i)][3]+shift,ANGS[1][str(i)][4]+shift,ANGS[1][str(i)][5]+shift))
            count+=1
    for mol in range(num_spec2+num_spec1, num_spec2+num_spec1+num_spec3):
        typstart = NTYPS[0][2] + NTYPS[1][2]
        for i in ANGS[2]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*num_spec2 + NCHAR[2][0]*(mol-num_spec1-num_spec2) +1
            dat.write("%s %s %s %s %s\n" % (count, typstart+ANGS[2][str(i)][0], ANGS[2][str(i)][3]+shift,ANGS[2][str(i)][4]+shift,ANGS[2][str(i)][5]+shift))
            count +=1 
    dat.write("\n")
    dat.close()

def write_dihedrals():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Dihedrals\n")
    dat.write("\n")
    count=1
    for mol in range(num_spec1):
        typstart = 0
        for i in DIHS[0]['name']:
            shift = NCHAR[0][0]*mol + 1
            dat.write("%s %s %s %s %s %s\n" % (count, typstart+DIHS[0][str(i)][0], DIHS[0][str(i)][4]+shift,DIHS[0][str(i)][5]+shift,DIHS[0][str(i)][6]+shift,DIHS[0][str(i)][7]+shift))
            count += 1
    for mol in range(num_spec2, num_spec2+num_spec1):
        typstart = NTYPS[0][3]
        for i in DIHS[1]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*(mol-num_spec1) + 1
            dat.write("%s %s %s %s %s %s\n" % (count, typstart+DIHS[1][str(i)][0], DIHS[1][str(i)][4]+shift,DIHS[1][str(i)][5]+shift,DIHS[1][str(i)][6]+shift,DIHS[0][str(i)][7]+shift))
            count += 1
    for mol in range(num_spec2+num_spec1, num_spec2+num_spec1+num_spec3):
        typstart = NTYPS[0][3] + NTYPS[1][3]
        for i in DIHS[2]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*num_spec2 + NCHAR[2][0]*(mol-num_spec1-num_spec2) +1
            dat.write("%s %s %s %s %s %s\n" % (count, typstart+DIHS[2][str(i)][0], DIHS[2][str(i)][4]+shift,DIHS[2][str(i)][5]+shift,DIHS[2][str(i)][6]+shift,DIHS[0][str(i)][7]+shift))
            count += 1
    dat.write("\n")
    dat.close()

def write_impropers():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Impropers\n")
    dat.write("\n")
    count=1
    for mol in range(num_spec1):
        typstart = 0
        for i in IMPS[0]['name']:
            shift = NCHAR[0][0]*mol + 1
            dat.write("%s %s %s %s %s %s\n" % (count, typstart+IMPS[0][str(i)][0], IMPS[0][str(i)][3]+shift,IMPS[0][str(i)][4]+shift,IMPS[0][str(i)][5]+shift,IMPS[0][str(i)][6]+shift))
            count += 1
    for mol in range(num_spec2, num_spec2+num_spec1):
        typstart = NTYPS[0][4]
        for i in IMPS[1]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*(mol-num_spec1) + 1
            dat.write("%s %s %s %s %s %s\n" % (count, typstart+DIHS[1][str(i)][0], IMPS[1][str(i)][3]+shift,IMPS[1][str(i)][4]+shift,IMPS[1][str(i)][5]+shift,IMPS[0][str(i)][6]+shift))
            count += 1
    for mol in range(num_spec2+num_spec1, num_spec2+num_spec1+num_spec3):
        typstart = NTYPS[0][4] + NTYPS[1][4]
        for i in IMPS[2]['name']:
            shift = NCHAR[0][0]*num_spec1 + NCHAR[1][0]*num_spec2 + NCHAR[2][0]*(mol-num_spec1-num_spec2) +1
            dat.write("%s %s %s %s %s %s\n" % (count, typstart+IMPS[2][str(i)][0], IMPS[2][str(i)][3]+shift,IMPS[2][str(i)][4]+shift,IMPS[2][str(i)][5]+shift,IMPS[0][str(i)][6]+shift))
            count += 1
    dat.write("\n") 
    dat.write("\n")
    dat.close()

def write_lammps_coeffs():
    lmpsfile = 'lmps.coeffs'
    lmps = open(lmpsfile, 'w')
    lmps.write("# Pair Coeffs species 1\n")
    a, indices = np.unique(ATMS[0]["atype"],return_index=True)
    count=1
    for j in indices:
        lmps.write("pair_coeff %s %s %9.3f %9.3f\n" % (ATMS[0]["atype"][j],ATMS[0]["atype"][j], ATMS[0]["kjeps"][j]/4.184, ATMS[0]["rmin"][j]/(2**(1/6.))))
        count += 1
    start = count-1
    lmps.write("# Pair Coeffs species 2\n")
    a, indices = np.unique(ATMS[1]["atype"],return_index=True)
    for j in indices:
        lmps.write("pair_coeff %s %s %9.3f %9.3f\n" % (ATMS[1]["atype"][j]+start,ATMS[1]["atype"][j]+start, ATMS[1]["kjeps"][j]/4.184, ATMS[1]["rmin"][j]/(2**(1/6.))))
        count += 1
    start = count-1
    lmps.write("# Pair Coeffs species 3\n")
    a, indices = np.unique(ATMS[2]["atype"],return_index=True)
    for j in indices:
        lmps.write("pair_coeff %s %s %9.3f %9.3f\n" % (ATMS[2]["atype"][j]+start,ATMS[2]["atype"][j]+start, ATMS[2]["kjeps"][j]/4.184, ATMS[2]["rmin"][j]/(2**(1/6.))))
        count += 1
    lmps.close()





NCHAR, NTYPS, ATMS, BNDS, ANGS, DIHS, IMPS = gen_packmol()
subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < IL_system.pmol"], shell=True)
x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
write_header()
write_atoms()
write_bonds()
write_angles()
write_dihedrals()
write_impropers()
write_lammps_coeffs()
