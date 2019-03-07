#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os, json

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
defaults        = {"nanion":0, "ncation":0, "nsolute":0, "tanion":"bmim", "tcation":"pf6", "tsolute":"R32","blength":40.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]),key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))
    
# Set parameters to input values or defaults
num_anion       = nml["nanion"]         if "nanion"         in nml else defaults["nanion"]
num_cation      = nml["ncation"]        if "ncation"        in nml else defaults["ncation"]
num_solute      = nml["nsolute"]        if "nsolute"        in nml else defaults["nsolute"]
typ_anion       = nml["tanion"]      if "tanion"      in nml else defaults["tanion"]
typ_cation      = nml["tcation"]     if "tcation"     in nml else defaults["tcation"]
typ_solute      = nml["tsolute"]     if "tsolute"     in nml else defaults["tsolute"]
blength         = nml["blength"]     if "blength"     in nml else defaults["blength"]

# Write out parameters
print("num_anion = %s of type %s" % (num_anion,typ_anion))
print("num_cation = %s of type %s" % (num_cation,typ_cation))
print("num_solute = %s of type %s" % (num_solute,typ_solute))
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
    natms=[]
    Q=[]
    MASS=[]
    TYP=[]
    if num_cation != 0:
        if typ_cation.lower() == "bmim":
            # This is the bmim structure in my lab notebook on p. 44
            bmim_file = "bmim.xyz"
            bmim = open(bmim_file, 'w')
            bmim.write("25\n")
            bmim.write("bmim\n")
            bmim.write("N1    0.000   -0.647   -0.704\n")
            bmim.write("N3   0.000    1.400    0.000\n")
            bmim.write("C2    0.000    0.624   -1.068\n")
            bmim.write("H1   -0.000    0.975   -2.079\n")
            bmim.write("C5    0.000   -0.681    0.676\n")
            bmim.write("H3    0.000   -1.601    1.219\n")
            bmim.write("C4    0.000    0.582    1.110\n")
            bmim.write("H2    0.000    0.974    2.104\n")
            bmim.write("C6   -0.000    2.872    0.038\n")
            bmim.write("H6   -0.879    3.217    0.567\n")
            bmim.write("H4    0.000    3.293   -0.957\n")
            bmim.write("H5    0.879    3.217    0.567\n")
            bmim.write("C7   -0.000   -1.861   -1.566\n")
            bmim.write("H8   -0.896   -2.421   -1.332\n")
            bmim.write("H7    0.896   -2.421   -1.332\n")
            bmim.write("C8    0.000   -1.591   -3.080\n")
            bmim.write("H9    0.876   -0.999   -3.339\n")
            bmim.write("H10   -0.876   -0.999   -3.339\n")
            bmim.write("C9    0.000   -2.870   -3.931\n")
            bmim.write("H11   -0.871   -3.477   -3.688\n")
            bmim.write("H12    0.871   -3.477   -3.688\n")
            bmim.write("C10    0.000   -2.579   -5.434\n")
            bmim.write("H13    0.877   -2.017   -5.736\n")
            bmim.write("H14   -0.8777  -2.017   -5.736\n")
            bmim.write("H15    0.000   -3.501   -6.007\n")
            bmim.close()
            q = [0.111, 0.133, 0.056, 0.177, -0.217, 0.207, -0.141, 0.181, -0.157, 0.142, 0.152, 0.073, 0.095, 0.045, 0.055, -0.122, 0.055, 0.001, 0.256, -0.029, -0.099, -0.209, 0.051, 0.040, 0.075]
            atype = [1, 1, 2, 3, 2, 4, 2, 4, 5, 6, 6, 6, 5, 7, 7, 8, 9, 9, 8, 9, 9, 10, 11, 11, 11]
            mass = [14.007, 12.011, 1.008, 1.008, 12.007, 1.008, 1.008, 12.011, 1.008, 12.011, 1.008]
            Q.append(q)
            MASS.append(mass)
            TYP.append(atype)
            # Write to the packmol file
            pm.write("structure bmim.xyz\n")
            pm.write("  number %s\n" % num_cation)
            pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
            pm.write("end structure\n")
            pm.write("\n")
            # Store the number of atoms
            natms.append(25)
            nbnds.append(25)
            nangs.append(35)
            ndihs.append(48)
            nimps.append(5)
            natyps.append(11)
            nbtyps.append(14)
            nangtyps.append(17)
            ndtyps.append(13)
            nityps.append(3)
        elif typ_cation.lower() == "emim":
            emim_file = "emim.xyz"
            emim = open(emim_file, 'w')
            emim.write("19\n")
            emim.write("emim\n")
            emim.write("C    0.890    0.908   -1.026\n")
            emim.write("C    2.944    1.764   -0.647\n")
            emim.write("C    2.070    2.771   -0.542\n")
            emim.write("H    0.088    0.234   -1.240\n")
            emim.write("H    4.006    1.818   -0.527\n")
            emim.write("H    2.290    3.794   -0.321\n")
            emim.write("C   -0.525    2.920   -0.782\n")
            emim.write("H   -1.312    2.262   -0.475\n")
            emim.write("H   -0.459    3.741   -0.100\n")
            emim.write("H   -0.733    3.288   -1.765\n")
            emim.write("N    0.749    2.187   -0.792\n")
            emim.write("C    2.688   -0.782   -1.168\n")
            emim.write("H    2.921   -0.919   -2.203\n")
            emim.write("C    3.966   -0.962   -0.327\n")
            emim.write("H    4.356   -1.947   -0.478\n")
            emim.write("H    4.696   -0.239   -0.628\n")
            emim.write("H    3.734   -0.825    0.709\n")
            emim.write("N    2.152    0.571   -0.961\n")
            emim.write("H    1.959   -1.504   -0.867\n")
            emim.close()
	    # Write to the packmol file
            pm.write("structure emim.xyz\n")
            pm.write("  number %s\n" % num_cation)
            pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
            pm.write("end structure\n")
            pm.write("\n")
            # Store the number of atoms
            natms.append(19)
        else:
            print("Chosen cation type is incorrect")
            sys.exit()
    if num_anion != 0:
        if typ_anion.lower() == "pf6":
            HFP_file = "PF6.xyz"
            HFP = open(HFP_file, 'w')
            HFP.write("7\n")
            HFP.write("PF6\n")
            HFP.write("P     0.000   0.000   0.000\n")
            HFP.write("F     0.000   1.646   0.000\n")
            HFP.write("F     0.000  -1.646   0.000\n")
            HFP.write("F     1.646   0.000   0.000\n")
            HFP.write("F    -1.646   0.000   0.000\n")
            HFP.write("F     0.000   0.000   1.646\n")
            HFP.write("F     0.000   0.000  -1.646\n")
            HFP.close()
            q = [ 1.458, -0.421, -0.426, -0.368, -0.368, -0.364, -0.414 ]
            atype = [1, 2, 2, 2, 2, 2, 2]
            mass  = [30.974, 18.997]
            Q.append(q)
            MASS.append(mass)
            TYP.append(atype)
            # Write to the packmol file
            pm.write("structure PF6.xyz\n")
            pm.write("  number %s\n" % num_anion)
            pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
            pm.write("end structure\n")
            pm.write("\n")
            natms.append(7)
            nbnds.append(6)
            nangs.append(12)
            ndihs.append(0)
            nimps.append(0)
            natyps.append(2)
            nbtyps.append(1)
            nangtyps.append(1)
            ndtyps.append(0)
            nityps.append(0)
        elif typ_anion.lower() == "bf4":
            BFF_file = "BF4.xyz"
            BFF = open(BFF_file, 'w')
            BFF.write("5\n")
            BFF.write("BF4\n")
            BFF.write("B    1.290   -0.392  -0.010\n")
            BFF.write("F    0.226   -1.449  -0.054\n")
            BFF.write("F    2.641   -1.042  -0.048\n")
            BFF.write("F    1.140   0.519   -1.191\n")
            BFF.write("F    1.152   0.402    1.255\n")
            BFF.close()
            q=[0.0,0.0,0.0,0.0,0.0]
            atype= [1,2,2,2,2]
            mass = [10.81, 18.998]
            MASS.append(mass)
            Q.append(q)
            TYP.append(atype)
            # Write to the packmol file
            pm.write("structure BF4.xyz\n")
            pm.write("  number %s\n" % num_anion)
            pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
            pm.write("end structure\n")
            pm.write("\n")
            natms.append(5)
            nbnds.append(4)
            natyps.append(2)
            nbtyps.append(1)
            nangtyps.append(1)
            ndtyps.append(0)
            nityps.append(0)
        else:
            print("Chosen anion type is incorrect")
            sys.exit()
    if num_solute != 0:
        if typ_solute.lower() == "r32":
            R32_file = "R32.xyz"
            R32 = open(R32_file, 'w')
            R32.write("5\n")
            R32.write("R32\n")
            R32.write("C     0.000    0.000   0.503\n")
            R32.write("F     0.000    1.110  -0.291\n")
            R32.write("F     0.000   -1.110  -0.291\n")
            R32.write("H    -0.908    0.000   1.107\n")
            R32.write("H     0.908    0.000   1.107\n")
            R32.close()
            atype = [1, 2, 2, 3, 3]
            q    = [ 0.0, 0.0, 0.0, 0.0, 0.0 ]
            mass = [12.011, 18.998, 1.008]
            Q.append(q)
            MASS.append(mass)
            TYP.append(atype)
            # Write to the packmol file
            pm.write("structure R32.xyz\n")
            pm.write("  number %s\n" % num_solute)
            pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
            pm.write("end structure\n")
            pm.write("\n")
            natms.append(5)
            nbnds.append(4)
            nangs.append(6)
            ndihs.append(0)
            nimps.append(0)
            natyps.append(3)
            nbtyps.append(2)
            nangtyps.append(3)
            ndtyps.append(0)
            nityps.append(0)
        else:
            print("Chosen solute type is incorrect")
            sys.exit()
    pm.close()
    return natms, Q, MASS, TYP

def write_header():
    datafile = "data.IL"
    dat      = open(datafile, 'w')
    dat.write("Ionic Liquid Data File for %s %s and %s\n" % (typ_anion, typ_cation, typ_solute))
    dat.write("\n")
    nATMS = natms[0]*num_cation + natms[1]*num_anion + natms[2]*num_solute
    nBNDS = nbnds[0]*num_cation + nbnds[1]*num_anion + nbnds[2]*num_solute
    nANGS = nangs[0]*num_cation + nangs[1]*num_anion + nangs[2]*num_solute
    nDIHS = ndihs[0]*num_cation + ndihs[1]*num_anion + ndihs[2]*num_solute
    nIMPS = nimps[0]*num_cation + nimps[1]*num_anion + nimps[2]*num_solute
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
    for i in range(len(MASS)):
        for j in range(len(MASS[i])):
            dat.write("%s %s\n" % (count, MASS[i][j]))
            count += 1
    dat.write("\n")

    dat.close()
    return

def write_atoms():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Atoms\n")
    dat.write("\n")
    for mol in range(num_cation+num_anion+num_solute):
        if mol < num_cation:
            moltype = 0
        elif mol >= num_cation and mol < num_cation+num_anion:
            moltype = 1
        else:
            moltype = 2
        for atom in range(natms[moltype]):
            mindex = mol + 1
            if moltype == 0:
                aindex   = mol*natms[moltype] + 1 + atom
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, TYP[moltype][atom], Q[moltype][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            elif moltype == 1:
                aindex = (mol-num_cation)*natms[moltype] + 1 + atom + num_cation*natms[0]
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, TYP[moltype][atom]+max(TYP[0]), Q[moltype][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
            else:
                aindex = (mol-num_cation-num_anion)*natms[moltype] + 1 + atom + num_cation*natms[0] + num_anion*natms[1]
                dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, TYP[moltype][atom]+max(TYP[0])+max(TYP[1]), Q[moltype][atom], x[aindex-1],y[aindex-1],z[aindex-1]))
    dat.write("\n")
    dat.close()
    return

def cation_bonding(atmstart, bndstart, angstart, dihstart, impstart, typstart, angtypstart, dihtypstart, imptypstart):
    line = [] # bond lines
    line2 = [] # angle lines
    line3 = [] # dihedral lines
    line4 = [] # improper lines

    if typ_cation.lower() == "bmim":
        """ This section is the bonding section and hopefully gets all the conectivity right 
            Bonding has been defined based on the paper Morrow and Maginn, J. Phys. Chem. B. 106, 2002."""
        N1,N3 =atmstart+1,atmstart+2
        C2,H1=atmstart+3,atmstart+4
        C5,H3=atmstart+5,atmstart+6
        C4,H2=atmstart+7,atmstart+8
        C6,H6,H4,H5=atmstart+9,atmstart+10,atmstart+11,atmstart+12
        C7,H8,H7=atmstart+13,atmstart+14,atmstart+15
        C8,H9,H10=atmstart+16,atmstart+17,atmstart+18
        C9,H11,H12=atmstart+19,atmstart+20,atmstart+21
        C10,H13,H14,H15=atmstart+22,atmstart+23,atmstart+24,atmstart+25
        line.append("%s %s %s %s\n" % (bndstart+1, typstart+1, N3, C6))         #N3-C6
        line.append("%s %s %s %s\n" % (bndstart+2, typstart+2, N1, C7))        #N1 C7
        line.append("%s %s %s %s\n" % (bndstart+3, typstart+3, C5, N1))         #C5-N1
        line.append("%s %s %s %s\n" % (bndstart+4, typstart+3, C4, N3))         #C4-C2
        line.append("%s %s %s %s\n" % (bndstart+5, typstart+4, C2, N1))         #C2-N1
        line.append("%s %s %s %s\n" % (bndstart+6, typstart+4, C2, N3))         #C2-N3
        line.append("%s %s %s %s\n" % (bndstart+7, typstart+5, C4, C5))         #C4-C5
        line.append("%s %s %s %s\n" % (bndstart+8, typstart+6, C2, H1))         #C2-H1
        line.append("%s %s %s %s\n" % (bndstart+9, typstart+7, H2, C4))         #H2-C4
        line.append("%s %s %s %s\n" % (bndstart+10, typstart+7, C5, H3))        #C5-H3
        line.append("%s %s %s %s\n" % (bndstart+11, typstart+8, C6, H6))       #C6-H6
        line.append("%s %s %s %s\n" % (bndstart+12, typstart+8, C6, H4))       #C6-H4
        line.append("%s %s %s %s\n" % (bndstart+13, typstart+8, C6, H5))       #C6-H5
        line.append("%s %s %s %s\n" % (bndstart+14, typstart+9, C7, H8))      #C7-H8
        line.append("%s %s %s %s\n" % (bndstart+15, typstart+9, C7, H7))      #C7-H7
        line.append("%s %s %s %s\n" % (bndstart+16, typstart+10, C8, H9))     #C8-H9
        line.append("%s %s %s %s\n" % (bndstart+17, typstart+10, C8, H10))     #C8-H10
        line.append("%s %s %s %s\n" % (bndstart+18, typstart+10, C9, H11))     #C9-H11
        line.append("%s %s %s %s\n" % (bndstart+19, typstart+10, C9, H12))     #C9-H12
        line.append("%s %s %s %s\n" % (bndstart+20, typstart+11, C10, H13))     #C10-H13
        line.append("%s %s %s %s\n" % (bndstart+21, typstart+11, C10, H14))     #C10-H14
        line.append("%s %s %s %s\n" % (bndstart+22, typstart+11, C10, H15))     #C10-H15
        line.append("%s %s %s %s\n" % (bndstart+23, typstart+12, C7, C8))     #C7-C8
        line.append("%s %s %s %s\n" % (bndstart+24, typstart+13, C8, C9))     #C8-C9
        line.append("%s %s %s %s\n" % (bndstart+25, typstart+14, C9, C10))     #C9-C10
        """ This section is the angle section and hopefully gets all the interior angles right. """
        line2.append("%s %s %s %s %s\n" % (angstart+1, angtypstart+1, C8, C7, N1)) 
        line2.append("%s %s %s %s %s\n" % (angstart+2, angtypstart+2, C5, N1, C2)) 
        line2.append("%s %s %s %s %s\n" % (angstart+3, angtypstart+2, C4, N3, C2)) 
        line2.append("%s %s %s %s %s\n" % (angstart+4, angtypstart+3, H4, C6, N3)) 
        line2.append("%s %s %s %s %s\n" % (angstart+5, angtypstart+3, H5, C6, N3)) 
        line2.append("%s %s %s %s %s\n" % (angstart+6, angtypstart+3, H6, C6, N3)) 
        line2.append("%s %s %s %s %s\n" % (angstart+7, angtypstart+4, H1, C2, N1))
        line2.append("%s %s %s %s %s\n" % (angstart+8, angtypstart+4, H1, C2, N3))
        line2.append("%s %s %s %s %s\n" % (angstart+9, angtypstart+5, N1, C5, C4))
        line2.append("%s %s %s %s %s\n" % (angstart+10, angtypstart+6, N1, C2, N3))
        line2.append("%s %s %s %s %s\n" % (angstart+11, angtypstart+7, N1, C2, H1))
        line2.append("%s %s %s %s %s\n" % (angstart+12, angtypstart+7, N3, C2, H1))
        line2.append("%s %s %s %s %s\n" % (angstart+13, angtypstart+8, H2, C4, C5))
        line2.append("%s %s %s %s %s\n" % (angstart+14, angtypstart+8, H3, C5, C4))
        line2.append("%s %s %s %s %s\n" % (angstart+15, angtypstart+9, N3, C4, H2))
        line2.append("%s %s %s %s %s\n" % (angstart+16, angtypstart+10, H4, C6, H5))
        line2.append("%s %s %s %s %s\n" % (angstart+17, angtypstart+10, H4, C6, H6))
        line2.append("%s %s %s %s %s\n" % (angstart+18, angtypstart+10, H5, C6, H5))
        line2.append("%s %s %s %s %s\n" % (angstart+19, angtypstart+11, H7, C7, H8))
        line2.append("%s %s %s %s %s\n" % (angstart+20, angtypstart+12, H7, C7, C8))
        line2.append("%s %s %s %s %s\n" % (angstart+21, angtypstart+12, H8, C7, C8))
        line2.append("%s %s %s %s %s\n" % (angstart+22, angtypstart+13, H9, C8, C7))
        line2.append("%s %s %s %s %s\n" % (angstart+23, angtypstart+13, H10, C8, C7))
        line2.append("%s %s %s %s %s\n" % (angstart+24, angtypstart+14, C7, C8, C9))
        line2.append("%s %s %s %s %s\n" % (angstart+25, angtypstart+14, C8, C9, C10))
        line2.append("%s %s %s %s %s\n" % (angstart+26, angtypstart+15, H9, C8, H10))
        line2.append("%s %s %s %s %s\n" % (angstart+27, angtypstart+15, H12, C9, H11))
        line2.append("%s %s %s %s %s\n" % (angstart+28, angtypstart+16, C8, C9, H11))
        line2.append("%s %s %s %s %s\n" % (angstart+29, angtypstart+16, C8, C9, H10))
        line2.append("%s %s %s %s %s\n" % (angstart+30, angtypstart+16, C9, C10, H13))
        line2.append("%s %s %s %s %s\n" % (angstart+31, angtypstart+16, C9, C10, H14))
        line2.append("%s %s %s %s %s\n" % (angstart+32, angtypstart+16, C9, C10, H15))
        line2.append("%s %s %s %s %s\n" % (angstart+33, angtypstart+17, H13, C10, H14))
        line2.append("%s %s %s %s %s\n" % (angstart+34, angtypstart+17, H13, C10, H15))
        line2.append("%s %s %s %s %s\n" % (angstart+35, angtypstart+17, H14, C10, H15))
        """ This section is the dihedral section and hopefully gets all the dihedrals right. """
        line3.append("%s %s %s %s %s %s\n" % (dihstart+1, dihtypstart+1, C2,N3,C4,C5))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+2, dihtypstart+1, N1,C5,C4,N3))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+3, dihtypstart+1, N1,C2,N3,C4))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+4, dihtypstart+2, H1,C2,N1,C5))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+5, dihtypstart+2, H1,C2,N3,C4))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+6, dihtypstart+3, H2,C4,C5,H3))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+7, dihtypstart+4, C4,C5,N1,C7))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+8, dihtypstart+4, C5,C4,N3,C6))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+9, dihtypstart+2, H2,C4,N3,C2))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+10, dihtypstart+2, N1,C5,C4,H2))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+11, dihtypstart+2, N3,C4,C5,H3))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+12, dihtypstart+5, N1,C2,N3,C6))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+13, dihtypstart+5, N3,C2,N1,C7))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+14, dihtypstart+5, H1,C2,N3,C6))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+15, dihtypstart+5, H1,C2,N1,C7))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+16, dihtypstart+5, H2,C4,N3,C6))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+17, dihtypstart+5, H3,C5,N1,C7))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+18, dihtypstart+6, C2,N1,C7,H8))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+19, dihtypstart+6, C2,N1,C7,H7))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+20, dihtypstart+6, C2,N3,C6,H4))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+21, dihtypstart+6, C2,N3,C6,H5))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+22, dihtypstart+6, C2,N3,C6,H6))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+23, dihtypstart+7, C4,N3,C6,H5))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+24, dihtypstart+7, C4,N3,C6,H4))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+25, dihtypstart+7, C4,N3,C6,H6))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+26, dihtypstart+7, C5,N1,C7,H7))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+27, dihtypstart+7, C5,N1,C7,H8))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+28, dihtypstart+8, C2,N1,C7,C8))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+29, dihtypstart+9, C5,N1,C7,C8))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+30, dihtypstart+10, N1,C7,C8,C9))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+31, dihtypstart+11, H12,C9,C10,H14))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+32, dihtypstart+11, H11,C9,C10,H14))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+33, dihtypstart+11, H12,C9,C10,H15))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+34, dihtypstart+11, H11,C9,C10,H15))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+35, dihtypstart+11, H12,C9,C10,H13))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+36, dihtypstart+11, H11,C9,C10,H13))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+37, dihtypstart+11, C8,C9,C10,H13))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+38, dihtypstart+11, C8,C9,C10,H14))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+39, dihtypstart+11, C8,C9,C10,H15))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+40, dihtypstart+10, N1,C7,C8,H9))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+41, dihtypstart+10, N1,C7,C8,H10))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+42, dihtypstart+12, C7,C8,C9,C10))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+43, dihtypstart+13, H7,C7,C8,H9))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+44, dihtypstart+13, H8,C7,C8,H10))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+45, dihtypstart+13, H10,C8,C9,H11))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+46, dihtypstart+13, H10,C8,C9,H12))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+47, dihtypstart+13, H9,C8,C9,H11))
        line3.append("%s %s %s %s %s %s\n" % (dihstart+48, dihtypstart+13, H9,C8,C9,H12))
        """ This section is the improper section and hopefully gets all the impropers right. """
        line4.append("%s %s %s %s %s %s\n" % (impstart+1, imptypstart+1, H1,C2,N3,N1))
        line4.append("%s %s %s %s %s %s\n" % (impstart+2, imptypstart+2, H2,C4,N3,C5))
        line4.append("%s %s %s %s %s %s\n" % (impstart+3, imptypstart+2, H3,C5,N1,C4))
        line4.append("%s %s %s %s %s %s\n" % (impstart+4, imptypstart+3, C7,N1,C2,C5))
        line4.append("%s %s %s %s %s %s\n" % (impstart+5, imptypstart+3, C6,N3,C2,C4))

    else: 
        print("Cation %s's bonding has not been defined yet... WHOOPS!" % typ_cation.lower())
    return line,line2,line3,line4
def anion_bonding(atmstart, bndstart, angstart, typstart,angtypstart):
    line = [] # bond lines
    line2 = [] # angle lines
    line3 = [] # dihedral lines
    line4 = [] # improper lines

    if typ_anion.lower() == "pf6":
        P1 = atmstart+1
        F1,F2,F3,F4,F5,F6 = atmstart+2, atmstart+3, atmstart+4, atmstart+5, atmstart+6, atmstart+7
        """ This section is the bonding section and hopefully gets all the conectivity right
            Bonding has been defined based on the paper Morrow and Maginn, J. Phys. Chem. B. 106, 2002."""
        line.append("%s %s %s %s\n" % (bndstart+1, typstart+1, P1, F1))
        line.append("%s %s %s %s\n" % (bndstart+2, typstart+1, P1, F2))
        line.append("%s %s %s %s\n" % (bndstart+3, typstart+1, P1, F3))
        line.append("%s %s %s %s\n" % (bndstart+4, typstart+1, P1, F4))
        line.append("%s %s %s %s\n" % (bndstart+5, typstart+1, P1, F5))
        line.append("%s %s %s %s\n" % (bndstart+6, typstart+1, P1, F6))
        """ This section is the angle section and hopefully gets all the interior angles right. """
        line2.append("%s %s %s %s %s\n" % (angstart+1, angtypstart+1, F1, P1, F2))
        line2.append("%s %s %s %s %s\n" % (angstart+2, angtypstart+1, F1, P1, F5))
        line2.append("%s %s %s %s %s\n" % (angstart+3, angtypstart+1, F1, P1, F4))
        line2.append("%s %s %s %s %s\n" % (angstart+4, angtypstart+1, F1, P1, F6))
        line2.append("%s %s %s %s %s\n" % (angstart+5, angtypstart+1, F2, P1, F5))
        line2.append("%s %s %s %s %s\n" % (angstart+6, angtypstart+1, F2, P1, F6))
        line2.append("%s %s %s %s %s\n" % (angstart+7, angtypstart+1, F5, P1, F4))
        line2.append("%s %s %s %s %s\n" % (angstart+8, angtypstart+1, F4, P1, F6))
        line2.append("%s %s %s %s %s\n" % (angstart+9, angtypstart+1, F3, P1, F5))
        line2.append("%s %s %s %s %s\n" % (angstart+10, angtypstart+1, F3, P1, F4))
        line2.append("%s %s %s %s %s\n" % (angstart+11, angtypstart+1, F3, P1, F6))
        line2.append("%s %s %s %s %s\n" % (angstart+12, angtypstart+1, F3, P1, F2))


    else:
        print("Anion %s 's bonding has not been defined yet... WHOOPS!" % typ_anion.lower())
    return line,line2

def solute_bonding(atmstart, bndstart, angstart, typstart, angtypstart):
    line = [] # bond lines
    line2 = [] # angle lines
    line3 = [] # dihedral lines
    line4 = [] # improper lines

    if typ_solute.lower() == "r32":
        C1=atmstart+1
        F1,F2 = atmstart+2, atmstart+3
        H1,H2 = atmstart+4, atmstart+5
        """ This section is the bonding section and hopefully gets all the conectivity right"""
        line.append("%s %s %s %s\n" % (bndstart+1, typstart+1, C1, F1))
        line.append("%s %s %s %s\n" % (bndstart+2, typstart+1, C1, F2))
        line.append("%s %s %s %s\n" % (bndstart+3, typstart+2, C1, H1))
        line.append("%s %s %s %s\n" % (bndstart+4, typstart+2, C1, H2))
        """ This section is the angle section and hopefully gets all the interior angles right. """
        line2.append("%s %s %s %s %s\n" % (angstart+1, angtypstart+1, F1, C1, F2))
        line2.append("%s %s %s %s %s\n" % (angstart+2, angtypstart+2, F1, C1, H1))
        line2.append("%s %s %s %s %s\n" % (angstart+3, angtypstart+2, F1, C1, H2))
        line2.append("%s %s %s %s %s\n" % (angstart+4, angtypstart+3, H1, C1, H2))
        line2.append("%s %s %s %s %s\n" % (angstart+5, angtypstart+2, H2, C1, F2))
        line2.append("%s %s %s %s %s\n" % (angstart+6, angtypstart+2, F2, C1, H1))
        
    else:
        print("Solute %s's bonding has not been defined yet... WHOOPS!" % typ_solute.lower())
    return line,line2


    
def write_bonds():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Bonds\n")
    dat.write("\n")
    for mol in range(num_cation):
        aindex   = mol*natms[0] + 1 
        bndstart = mol*nbnds[0]
        typstart = 0
        line, line2, line3, line4 = cation_bonding(aindex-1, bndstart,0,0,0, typstart,0,0,0)
        for i in range(nbnds[0]):
            dat.write(line[i])
    for mol in range(num_cation, num_cation+num_anion):
        aindex   = (mol-num_cation)*natms[1] + 1 + num_cation*natms[0]
        bndstart = (mol-num_cation)*nbnds[1] + num_cation*nbnds[0]
        typstart = nbtyps[0]
        line, line2= anion_bonding(aindex-1, bndstart,0, typstart,0)
        for i in range(nbnds[1]):
            dat.write(line[i])
    for mol in range(num_cation+num_anion, num_cation+num_anion+num_solute):
        aindex   = (mol-num_cation-num_anion)*natms[2] + 1 + num_cation*natms[0] + num_anion*natms[1]
        bndstart = (mol-num_cation-num_anion)*nbnds[2] + num_cation*nbnds[0] + num_anion*nbnds[1]
        typstart = nbtyps[0] + nbtyps[1]
        line, line2 = solute_bonding(aindex-1, bndstart,0, typstart,0)
        for i in range(nbnds[2]):
            dat.write(line[i])
    dat.write("\n")
    dat.close()

def write_angles():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Angles\n")
    dat.write("\n")
    for mol in range(num_cation):
        aindex = mol*natms[0] + 1
        bndstart, angstart = mol*nbnds[0], mol*nangs[0]
        typstart = 0
        angtypstart = 0
        line, line2,line3,line4 = cation_bonding(aindex-1, bndstart, angstart,0,0, typstart, angtypstart,0,0)
        for i in range(nangs[0]):
            dat.write(line2[i])
    for mol in range(num_cation, num_cation+num_anion):
        aindex   = (mol-num_cation)*natms[1] + 1 + num_cation*natms[0]
        bndstart, angstart = (mol-num_cation)*nbnds[1] + num_cation*nbnds[0], (mol-num_cation)*nangs[1] + num_cation*nangs[0]
        typstart = nbtyps[0]
        angtypstart = nangtyps[0]
        print(angtypstart)
        line,line2 = anion_bonding(aindex-1, bndstart, angstart, typstart, angtypstart)
        for i in range(nangs[1]):
            dat.write(line2[i])
    for mol in range(num_cation+num_anion, num_cation+num_anion+num_solute):
        aindex   = (mol-num_cation-num_anion)*natms[2] + 1 + num_cation*natms[0] + num_anion*natms[1]
        bndstart,angstart = (mol-num_cation-num_anion)*nbnds[2] + num_cation*nbnds[0] + num_anion*nbnds[1],(mol-num_cation-num_anion)*nangs[2] + num_cation*nangs[0] + num_anion*nangs[1]
        typstart = nbtyps[0] + nbtyps[1]
        angtypstart = nangtyps[0]+nangtyps[1]
        line, line2 = solute_bonding(aindex-1, bndstart,angstart, typstart,angtypstart)
        for i in range(nangs[2]):
            dat.write(line2[i])
    dat.write("\n")
    dat.close()

def write_dihedrals():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Dihedrals\n")
    dat.write("\n")
    for mol in range(num_cation):
        aindex = mol*natms[0] + 1
        bndstart, angstart, dihstart= mol*nbnds[0], mol*nangs[0], mol*ndihs[0]
        typstart = 0
        angtypstart = 0
        dihtypstart = 0
        line, line2, line3, line4 = cation_bonding(aindex-1, bndstart, angstart,dihstart,0, typstart, angtypstart,dihtypstart,0)
        for i in range(ndihs[0]):
            dat.write(line3[i])
    dat.write("\n")
    dat.close()

def write_impropers():
    datafile = "data.IL"
    dat      = open(datafile, 'a')
    dat.write("Impropers\n")
    dat.write("\n")
    for mol in range(num_cation):
        aindex = mol*natms[0] + 1
        bndstart, angstart, dihstart, impstart = mol*nbnds[0], mol*nangs[0], mol*ndihs[0], mol*nimps[0]
        typstart = 0
        angtypstart = 0
        dihtypstart = 0
        imptypstart = 0
        line, line2, line3, line4 = cation_bonding(aindex-1, bndstart, angstart,dihstart,impstart, typstart, angtypstart,dihtypstart,imptypstart)
        for i in range(nimps[0]):
            dat.write(line4[i])
    dat.write("\n")
    dat.close()



# Initialize Arrays
natyps=[]
nbtyps=[]
nangtyps=[]
ndtyps=[]
nityps=[]
natms=[]
nbnds=[]
nangs=[]
ndihs=[]
nimps=[]

natms, Q , MASS, TYP = gen_packmol()
print(natms)
subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < IL_system.pmol"], shell=True)
x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
write_header()
write_atoms()
write_bonds()
write_angles()
write_dihedrals()
write_impropers()
