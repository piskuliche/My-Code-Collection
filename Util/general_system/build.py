#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os, json
import molecular_species as molspec

"""
This is a python program which can build a box of generally defined molecules.
These molecules cane be included in molecular_species.py


NOTE:
    Pair-Styles:
        lj rmin form or lj 4eps form
    Bond-Styles: 
        Rigid (RASPA)
        Harmonic
    Angle Style 
        Rigid (RASPA)
        Harmonic 
    Dihedral Style 
        Rigid 
        CVFF (K*(1+cos(n*x)))
        Charmm (K*(1+cos(n*x-theta)))
    Improper Style
        Rigid
        CVFF 



"""

if not os.path.exists('./tmp'):
    os.makedirs('./tmp')

if not os.path.exists('./raspa'):
    os.makedirs('./raspa')

if not os.path.exists('./connect'):
    os.makedirs('./connect')

# Read the NML file
try:
    nml = json.load(sys.stdin)
except ValueError:
#json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys, and typcheck values
defaults        = {"nspec":[0], "tspec":["bmim"],"blength":40.0,"num_components":1, "ff_type":["lj"], "e_unit":["kcal"],"f_unit":["kcal"], "eo_unit":["kcal"],"fo_unit":["kcal"], "shift_f":1.0}
    
# Set parameters to input values or defaults
num_spec       = nml["nspec"]          if "nspec"          in nml else defaults["nspec"]
typ_spec       = nml["tspec"]          if "tspec"          in nml else defaults["tspec"]
blength         = nml["blength"]        if "blength"        in nml else defaults["blength"]
num_components  = nml["num_components"] if "num_components" in nml else defaults["num_components"]
fftype          = nml["ff_type"]       if "ff_type"       in nml else defaults["ff_type"]
eunit           = nml["e_unit"]         if "e_unit"         in nml else defaults["e_unit"]
funit           = nml["f_unit"]         if "f_unit"         in nml else defaults["f_unit"]
eounit          = nml["eo_unit"]         if "eo_unit"         in nml else defaults["eo_unit"]
founit          = nml["fo_unit"]         if "fo_unit"         in nml else defaults["fo_unit"]
fshift          = nml["shift_f"]         if "shift_f"         in nml else defaults["shift_f"]

# Write out parameters
print("There are %s total components" % num_components)
for i in range(num_components):
    print("Species %s has %s %s molecules." % (i, num_spec[i],typ_spec[i]))
print("box length = %s" % blength)


# Conversions
conv_to_kcal = { 'kj':0.2390057, 'K':0.0019872041, 'kcal':1.0, 'ev':627.509}
conv_energy=[]
conv_force=[]
for i in range(num_components):
    conv_energy.append(conv_to_kcal[eunit[i]] * 1.0/conv_to_kcal[eounit[i]])
    conv_force.append(conv_to_kcal[funit[i]]  * 1.0/conv_to_kcal[founit[i]] * fshift)

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
                dat.write("%s %s %s %s %s %s\n" % (count, typstart+IMPS[i][str(imp)][0], IMPS[i][str(imp)][4]+shiftimps,IMPS[i][str(imp)][5]+shiftimps,IMPS[i][str(imp)][6]+shiftimps,IMPS[i][str(imp)][7]+shiftimps))
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
            if fftype[i] == "rmin":
                    lmps.write("pair_coeff %s %s %9.5f %9.5f\n" % (ATMS[i]["atype"][j]+start,ATMS[i]["atype"][j]+start, ATMS[i]["eps"][j]*conv_energy[i], ATMS[i]["rmin"][j]/(2**(1/6.))))
            elif fftype[i] == "lj":
                    lmps.write("pair_coeff %s %s %9.5f %9.5f\n" % (ATMS[i]["atype"][j]+start,ATMS[i]["atype"][j]+start, ATMS[i]["eps"][j]*conv_energy[i], ATMS[i]["sig"][j]))
            else:
                print("Unrecognized force field type.")
                
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
                lmps.write("bond_coeff %s %9.3f %9.3f\n" % (BNDS[i][bond][0]+start, BNDS[i][bond][1]*conv_force[i], BNDS[i][bond][2]))
                usedtyps.append(BNDS[i][bond][0]+start)
        if usedtyps:
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
                lmps.write("angle_coeff %s %9.5f %9.5f\n" % (ANGS[i][angle][0]+start, ANGS[i][angle][1]*conv_force[i], ANGS[i][angle][2]))
                usedtyps.append(ANGS[i][angle][0]+start)
        if usedtyps:
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
                lmps.write("dihedral_coeff %s %9.5f %s %d 0.0\n" % (DIHS[i][dihedral][0]+start, DIHS[i][dihedral][1]*conv_force[i], DIHS[i][dihedral][2], DIHS[i][dihedral][3]))
                usedtyps.append(DIHS[i][dihedral][0]+start)
        if usedtyps:
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
                lmps.write("improper_coeff %s %9.5f %9.5f\n" % (IMPS[i][improper][0]+start, IMPS[i][improper][1]*conv_force[i], IMPS[i][improper][2]))
                usedtyps.append(IMPS[i][improper][0]+start)
        if usedtyps:
            start = max(usedtyps)
    lmps.close()
    return


def write_raspa_ff():
    # write the pseudo atom file
    pseudo_file = "./raspa/pseudo_atoms.def"
    pafile = open(pseudo_file, 'w')
    count = 0
    for i in range(num_components):
        for name in np.unique(ATMS[i]['name']):
            count += 1
    pafile.write("#number of pseudo atoms\n")
    pafile.write("%s\n" % count)
    pafile.write("#type print as scat oxidation mass charge polarization B-factor radii connectivity anisotropic anisotropic-type tinker-type\n")
    for i in range(num_components):
        for name in np.unique(ATMS[i]['name']):
            index = ATMS[i]['name'].index(name)
            pafile.write("%s_%s yes %s  %s  0 %9.6f %9.6f 0.0 0.0 0.0 0 0.0 absolute 0\n" % (typ_spec[i].lower(), name, name[:1], name[:1], ATMS[i]['mass'][index], ATMS[i]['q'][index]))
    pafile.close()

    # Write the ff definition file
    def_file = "./raspa/force_field.def"
    dfile = open(def_file, 'w')
    dfile.write("# rules to overwrite\n")
    dfile.write("0\n")
    dfile.write("# of defined interactions\n")
    dfile.write("# mixing rules to overwrite\n")
    dfile.write("0\n")
    dfile.close()
    # Write the ff mixing rules file
    mr_raspa_file = "./raspa/force_field_mixing_rules.def"
    mrfile = open(mr_raspa_file, 'w')
    mrfile.write("# general rule for shifted vs truncated\n")
    mrfile.write("truncated\n")
    mrfile.write("#general rule tailcorrections\n")
    mrfile.write("yes\n")
    mrfile.write("#number of defined interactions\n")
    count = 0
    for i in range(num_components):
        for name in np.unique(ATMS[i]['name']):
            count += 1
    mrfile.write("%s\n" % count)
    mrfile.write("# type interaction\n")
    for i in range(num_components):
        for name in np.unique(ATMS[i]['name']):
            index = ATMS[i]['name'].index(name)
            if fftype[i] == "rmin":
                mrfile.write("%s_%s LENNARD_JONES %9.4f %9.4f\n" % (typ_spec[i].lower(), name, ATMS[i]['eps'][index]*conv_energy[i],  ATMS[i]['rmin'][index]/(2**(1/6.))))
            elif fftype[i] == "lj":
                mrfile.write("%s_%s LENNARD_JONES %9.4f %9.4f\n" % (typ_spec[i].lower(), name, ATMS[i]['eps'][index]*conv_energy[i],  ATMS[i]['sig'][index]/(2**(1/6.))))

    mrfile.write("# general mixing rule for Lennard-Jones\n")
    mrfile.write("Lorentz-Berthelot\n")
    mrfile.close()
    # Write the molecule definition file
    vdwpairs = Intra_VDW()
    shiftatoms = 0
    for i in range(num_components):
        molfile = './raspa/'+typ_spec[i].lower()+'.def'
        mol = open(molfile, 'w')
        mol.write('# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]\n')
        mol.write('300.0\n')
        mol.write('4000000\n')
        mol.write('0.393\n')
        mol.write('# Number of Atoms\n')
        mol.write('%s\n' % NCHAR[i][0])
        mol.write('# of Groups\n')
        mol.write('%s\n' % len(np.unique(ATMS[i]['group'])))
        count = 0
        for j in np.unique(ATMS[i]['group']):
            # Rigid
            if ATMS[i]['typgrp'][ATMS[i]['group'].index(j)] == 0:
                mol.write('# Rigid Group\n')
                mol.write('rigid\n')
                mol.write('# Number of atoms\n')
                count2=0
                for atom in ATMS[i]['name']:
                    if ATMS[i]['group'][ATMS[i]['name'].index(atom)] == j:
                        count2 += 1
                mol.write('%s\n' % count2)
                mol.write('# atomic positions\n')
                for atom in ATMS[i]['name']:
                    if ATMS[i]['group'][ATMS[i]['name'].index(atom)] == j:
                        aindex = 1 + ATMS[i]['name'].index(atom) + shiftatoms
                        mol.write('%s %s_%s %9.5f %9.5f %9.5f\n' % (ATMS[i]['name'].index(atom), typ_spec[i], atom, x[aindex-1], y[aindex-1], z[aindex-1]))
                        count += 1
            # Flexible
            if ATMS[i]['typgrp'][ATMS[i]['group'].index(j)] == 1:
                mol.write('# Flexible Group\n')
                mol.write('flexible\n')
                mol.write('# Number of atoms\n')
                count2=0
                for atom in ATMS[i]['name']:
                    if ATMS[i]['group'][ATMS[i]['name'].index(atom)] == j:
                        count2 += 1
                mol.write('%s\n' % count2)
                mol.write('# atomic positions\n')
                for atom in ATMS[i]['name']:
                    if ATMS[i]['group'][ATMS[i]['name'].index(atom)] == j:
                        mol.write('%s %s_%s\n' % (ATMS[i]['name'].index(atom), typ_spec[i], atom))
                        count += 1
        shiftatoms += num_spec[i]*NCHAR[i][0]
        mol.write('# Chiral centers Bond Bond Dipoles Bend UrayBradley InvBend Torsion Imp. Torsion Bond/Bond Stretch/Bend Bend/Bend Stretch/Torsion Bend/Torsion IntraVDW IntraCoulomb\n')
        mol.write("0 %s 0 %s 0 0 %s %s 0 0 0 0 0 %s 0\n" % (NCHAR[i][1],NCHAR[i][2],NCHAR[i][3],NCHAR[i][4],len(vdwpairs[i])))
        mol.write("Bond Stretch: atom n1-n2, type, parameters\n")
        for bond in BNDS[i]['name']:
            if ATMS[i]['typgrp'][BNDS[i][bond][3]] == 0 and ATMS[i]['typgrp'][BNDS[i][bond][4]] == 0:
                mol.write("%s %s RIGID_BOND\n" % (BNDS[i][bond][3],BNDS[i][bond][4]))
            else:
                mol.write("%s %s HARMONIC_BOND %9.3f %9.3f\n" % (BNDS[i][bond][3],BNDS[i][bond][4],BNDS[i][bond][1]*conv_force[i],BNDS[i][bond][2]))
        mol.write("# Angle bending: atom n1-n2-n3, type, parameters\n")
        for angle in ANGS[i]['name']:
            if ATMS[i]['typgrp'][ANGS[i][angle][3]] == 0 and ATMS[i]['typgrp'][ANGS[i][angle][4]] == 0 and ATMS[i]['typgrp'][ANGS[i][angle][5]] == 0:
                mol.write("%s %s %s RIGID_BEND\n" % (ANGS[i][angle][3],ANGS[i][angle][4], ANGS[i][angle][5]))
            else:
                mol.write("%s %s %s HARMONIC_BEND %9.3f %9.3f\n" % (ANGS[i][angle][3],ANGS[i][angle][4],ANGS[i][angle][5], ANGS[i][angle][1]*conv_force[i],ANGS[i][angle][2]))
        mol.write("# Torsion n1-n2-n3-n4 type\n")
        for dihedral in DIHS[i]['name']:
            if ATMS[i]['typgrp'][DIHS[i][dihedral][4]] == 0 and ATMS[i]['typgrp'][DIHS[i][dihedral][5]] == 0 and ATMS[i]['typgrp'][DIHS[i][dihedral][6]] == 0 and ATMS[i]['typgrp'][DIHS[i][dihedral][7]] == 0:
                mol.write("%s %s %s %s RIGID_DIHEDRAL\n" % (DIHS[i][dihedral][4],DIHS[i][dihedral][5], DIHS[i][dihedral][6],DIHS[i][dihedral][7]))
            else:
                mol.write("%s %s %s %s CVFF_DIHEDRAL %9.3f %s %9.3f\n" % (DIHS[i][dihedral][4],DIHS[i][dihedral][5],DIHS[i][dihedral][6],DIHS[i][dihedral][7],DIHS[i][dihedral][1]*conv_force[i],DIHS[i][dihedral][2], DIHS[i][dihedral][3]))
        mol.write("# Improper Torsion n1-n2-n3-n4 type\n")
        for improper in IMPS[i]['name']:
            if ATMS[i]['typgrp'][IMPS[i][improper][4]] == 0 and ATMS[i]['typgrp'][IMPS[i][improper][5]] == 0 and ATMS[i]['typgrp'][IMPS[i][improper][6]] == 0 and ATMS[i]['typgrp'][IMPS[i][improper][7]] == 0:
                mol.write("%s %s %s %s RIGID_IMPROPER\n" % (IMPS[i][improper][4],IMPS[i][improper][5], IMPS[i][improper][6], IMPS[i][improper][7]))
            else:
                mol.write("%s %s %s %s CVFF_IMPROPER %9.3f %s %9.3f\n" % (IMPS[i][improper][4],IMPS[i][improper][5],IMPS[i][improper][6],IMPS[i][improper][7], IMPS[i][improper][1]*conv_force[i],IMPS[i][improper][2], IMPS[i][improper][3]))
        mol.write("# Intra VDW: atom n1-n2\n")
        for pair in vdwpairs[i]:
            mol.write("%s %s\n" % (pair[0],pair[1]))
        mol.write("# Number of config moves\n")
        mol.write("%s\n" % len(ATMS[i]['cf']))
        mol.write("# nr fixed, list\n")
        for item in ATMS[i]['cf']:
            mol.write("%s\n" % item)
        mol.close()
    return

def Intra_VDW():
    InVDW=[]
    VDWALLPAIRS=[]
    for i in range(num_components):
        VDWPAIR = []
        for atom in ATMS[i]['name']:
            index = ATMS[i]['name'].index(atom)
            line = []
            if len(DIHS[i]['name']) != 0:
                for dihedral in DIHS[i]['name']:
                    if index == DIHS[i][dihedral][4]:
                        line.append(DIHS[i][dihedral][5])
                        line.append(DIHS[i][dihedral][6])
                        line.append(DIHS[i][dihedral][7])
                    elif index == DIHS[i][dihedral][5]:
                        line.append(DIHS[i][dihedral][4])
                        line.append(DIHS[i][dihedral][6])
                        line.append(DIHS[i][dihedral][7]) 
                    elif index == DIHS[i][dihedral][6]:
                        line.append(DIHS[i][dihedral][4])
                        line.append(DIHS[i][dihedral][5])
                        line.append(DIHS[i][dihedral][7])
                    elif index == DIHS[i][dihedral][7]:
                        line.append(DIHS[i][dihedral][4])
                        line.append(DIHS[i][dihedral][5])
                        line.append(DIHS[i][dihedral][6])
                for atm in range(len(ATMS[i]['name'])):
                    if atm not in line and [atm, index] not in VDWPAIR and atm != index:
                        VDWPAIR.append([index, atm])
        InVDW.append(len(VDWPAIR))
        print("There are %s IntraVDW interactions" % InVDW[i])
        VDWALLPAIRS.append(VDWPAIR)
    return VDWALLPAIRS
            

def write_connectivity():
    """ This writes a barebone connectivity file readable by lammps as described on
        https://lammps.sandia.gov/doc/molecule.html for the molecule command. Instead
        of coordinates it leaves it blank."""
    typ_start =[0,0,0,0,0]
    for i in range(num_components):
        connfile = "./connect/"+typ_spec[i].lower()+".connect"
        cnf = open(connfile,'w')
        cnf.write("# This file is a connectivity file for %s\n" % typ_spec[i].lower())
        # Header Section
        cnf.write("%d atoms\n" % NCHAR[i][0])
        cnf.write("%d bonds\n" % NCHAR[i][1])
        cnf.write("%d angles\n"% NCHAR[i][2])
        cnf.write("%d dihedrals\n"% NCHAR[i][3])
        cnf.write("%d impropers\n"% NCHAR[i][4])
        cnf.write("\n")
        # Coordinates
        cnf.write("Coords\n")
        cnf.write("\n")
        for atm in range(NCHAR[i][0]):
            cnf.write("%d\n" % (atm+1))
        cnf.write("\n")
        # Types
        cnf.write("Types\n")
        cnf.write("\n")
        for atm in range(NCHAR[i][0]):
            cnf.write("%d %d\n" % (atm+1,ATMS[i]["atype"][atm]+typ_start[0]))
        cnf.write("\n")
        cnf.write("Charges\n")
        cnf.write("\n")
        for atm in range(NCHAR[i][0]): 
            cnf.write("%d %.4f\n" % (atm+1,ATMS[i]["q"][atm])) 
        cnf.write("\n")
        cnf.write("Masses\n")
        cnf.write("\n")
        for atm in range(NCHAR[i][0]):        
            cnf.write("%d %.4f\n" % (atm+1,ATMS[i]["mass"][atm]))                                   
        cnf.write("\n")
        if NCHAR[i][1] != 0:
            cnf.write("Bonds\n")
            cnf.write("\n")
            count = 1
            for bnd in BNDS[i]["name"]:        
                cnf.write("%d %d %d %d\n" % (count,BNDS[i][str(bnd)][0]+typ_start[1],BNDS[i][str(bnd)][3]+1,BNDS[i][str(bnd)][4]+1))
                count += 1                                  
            cnf.write("\n")
        if NCHAR[i][2] != 0:
            cnf.write("Angles\n")
            cnf.write("\n") 
            count = 1
            for ang in ANGS[i]["name"]:      
                cnf.write("%d %d %d %d %d\n" % (count,ANGS[i][str(ang)][0]+typ_start[2],ANGS[i][str(ang)][3]+1,ANGS[i][str(ang)][4]+1,ANGS[i][str(ang)][5]+1))                                                   
                count += 1
            cnf.write("\n")
        if NCHAR[i][3] != 0:
            cnf.write("Dihedrals\n")
            cnf.write("\n") 
            count = 1
            for dih in DIHS[i]["name"]:      
                cnf.write("%d %d %d %d %d %d\n" % (count,DIHS[i][str(dih)][0]+typ_start[3],DIHS[i][str(dih)][4]+1,DIHS[i][str(dih)][5]+1,DIHS[i][str(dih)][6]+1,DIHS[i][str(dih)][7]+1))    
                count += 1
            cnf.write("\n")
        if NCHAR[i][4] != 0:
            cnf.write("Impropers\n")
            cnf.write("\n") 
            count = 1
            for imp in IMPS[i]["name"]:      
                cnf.write("%d %d %d %d %d %d\n" % (count,IMPS[i][str(imp)][0]+typ_start[4],IMPS[i][str(imp)][4]+1,IMPS[i][str(imp)][5]+1,IMPS[i][str(imp)][6]+1,IMPS[i][str(imp)][7]+1))          
                count += 1
        cnf.write("\n")


        for typ in range(5):
            typ_start[typ] += NTYPS[i][typ]
        cnf.close()
    return


def connectflag(NCHAR):
    """Small function that checks if there are each type of connectivity"""
    flags = [False,False,False,False,False]
    for i in range(5):
        if np.sum(np.transpose(NCHAR)[i]) != 0:
            flags[i]=True
    return flags




NCHAR, NTYPS, ATMS, BNDS, ANGS, DIHS, IMPS = gen_packmol()
subprocess.call(["/panfs/pfs.local/home/s393s653/programs/Executable/packmol < system.pmol > packmol.output"], shell=True)
x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
flags = connectflag(NCHAR)
write_header()
if flags[0] == True:
    write_atoms()
if flags[1] == True:
    write_bonds()
if flags[2] == True:
    write_angles()
if flags[3] == True:
    write_dihedrals()
if flags[4] == True:
    write_impropers()
write_lammps_paircoeffs()
write_lammps_bondcoeffs()
write_lammps_anglecoeffs()
write_lammps_dihedralcoeffs()
write_lammps_impropercoeffs()
#write_raspa_ff()
write_connectivity()
