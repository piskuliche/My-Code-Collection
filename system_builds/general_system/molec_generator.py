#!/usr/bin/env python

def read_connect_file(filename):
    """ This is a subroutine to read a lammps-style connectivity file as described on https://lammps.sandia.gov/doc/molecule.html"""
    natms, nbnds, nangs, ndihs, nimps = 0,0,0,0,0
    fCoords,    cCoords     = 0, 0
    fTypes,     cTypes      = 0, 0
    fCharges,   cCharges    = 0, 0
    fMasses,    cMasses     = 0, 0
    fBonds,     cBonds      = 0, 0
    fAngles,    cAngles     = 0, 0
    fDihedrals, cDihedrals  = 0, 0 
    fImpropers, cImpropers  = 0, 0
    coords      = []
    types       = []
    charges     = []
    masses      = []
    bonds       = []
    angles      = []
    dihedrals   = []
    impropers   = []
    
    print("Reading connectivity file")
    with open(filename) as f:
        for line in f:
            if "atoms" in line:
                natms = int(line.split()[0])
            elif "bonds" in line:
                nbnds = int(line.split()[0])
            elif "angles" in line:
                nangs = int(line.split()[0])
            elif "dihedrals" in line:
                ndihs = int(line.split()[0])
            elif "impropers" in line:
                nimps = int(line.split()[0])
            if "Coords" in line:
                fCoords = 1
            if fCoords == 1:
                cCoords += 1
            if fCoords == 1 and cCoords > 2 and cCoords <= 2 + natms:
                coords.append([int(line.split()[0]),float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            elif cCoords > 2 + natms:
                fCoords = 0
            if "Types" in line:
                fTypes = 1
            if fTypes == 1:
                cTypes += 1
            if fTypes == 1 and cTypes > 2 and cTypes <= 2 + natms:
                types.append([int(line.split()[0]), int(line.split()[1])])
            elif cTypes > 2 + natms:
                fTypes = 0
            if "Charges" in line:
                fCharges = 1
            if fCharges == 1:
                cCharges += 1
            if fCharges == 1 and cCharges > 2 and cCharges <= 2 + natms:
                charges.append([int(line.split()[0]), float(line.split()[1])])
            elif cCharges > 2 + natms:  
                fCharges = 0
            if "Masses" in line:
                fMasses = 1
            if fMasses == 1:
                cMasses += 1
            if fMasses == 1 and cMasses > 2 and cMasses <= 2 + natms:
                masses.append([int(line.split()[0]), float(line.split()[1])])
            elif cMasses > 2 + natms:
                fMasses = 0
            if "Bonds" in line:
                fBonds = 1
            if fBonds == 1:
                cBonds += 1
            if fBonds == 1 and cBonds > 2 and cBonds <= 2 + nbnds:
                bonds.append([int(line.split()[0]), int(line.split()[1]),int(line.split()[2]), int(line.split()[3])])
            elif cBonds > 2 + nbnds:
                fBonds = 0
            if "Angles" in line:
                fAngles = 1
            if fAngles == 1:
                cAngles += 1
            if fAngles == 1 and cAngles > 2 and cAngles <= 2 + nangs:
                angles.append([int(line.split()[0]), int(line.split()[1]),int(line.split()[2]), int(line.split()[3]), int(line.split()[4])])
            elif cAngles > 2 + nangs:
                fAngles = 0
            if "Dihedrals" in line:
                fDihedrals = 1
            if fDihedrals == 1:
                cDihedrals += 1
            if fDihedrals == 1 and cDihedrals > 2 and cDihedrals <= 2 + ndihs:
                dihedrals.append([int(line.split()[0]), int(line.split()[1]),int(line.split()[2]), int(line.split()[3]), int(line.split()[4]), int(line.split()[5])])
            elif cDihedrals > 2 + ndihs:
                fDihedrals = 0
            if "Impropers" in  line:
                fImpropers = 1
            if fImpropers == 1:
                cImpropers += 1
            if fImpropers == 1 and cImpropers > 2 and cImpropers <= 2 + nimps:
                impropers.append([int(line.split()[0]), int(line.split()[1]),int(line.split()[2]), int(line.split()[3]), int(line.split()[4]), int(line.split()[5])])
            elif cImpropers > 2 + nimps:
                fImpropers = 0
    nchar = [natms, nbnds, nangs, ndihs, nimps]
    print("Read Connectivity file")
    return nchar, coords, types, charges, masses, bonds, angles, dihedrals, impropers

def calc_ntypes(types, bonds, angles, dihedrals, impropers):
    ntyps = []
    ntyps.append(np.max(np.transpose(types)[1]))
    if len(bonds) > 0:
        ntyps.append(np.max(np.transpose(bonds)[1]))
    if len(angles) > 0:
        ntyps.append(np.max(np.transpose(angles)[1]))
    if len(dihedrals) > 0:
        ntyps.append(np.max(np.transpose(dihedrals)[1]))
    if len(impropers) > 0:
        ntyps.append(np.max(np.transpose(impropers)[1]))
    return ntyps

def write_atm_line(secname, array,nchar):
    line = ("        '%s':[" % secname)
    for i in range(nchar[0]):
        if i != 0:
            line += (',"%s"' % array[i])
        else:
            line += ('"%s"' % array[i])
    line += ("],\n")
    return line

def read_names(molname):
    namefile = str(molname)+".names"
    names = []
    with open(namefile) as nam:
        for line in nam:
            names.append(str(line.split()[1]))
    return names

def read_paircoeffs():
    pcoeffs = []
    with open('lmps.paircoeffs') as lmps:
        for line in lmps:
            pcoeffs.append(int(line.split()[1]), int(line.split()[2]), float(line.split(3)), float(line.split(4)))
    return pcoeffs

def write_molec_py(infile,outfile,molname):
    nchar, coords, types, charges, masses, bonds, angles, dihedrals, impropers = read_connect_file(infile)
    ntyps = calc_ntypes(types, bonds, angles, dihedrals, impropers)
    names = read_names(molname)
    zeros = np.zeros(nchar[0])

    py = open(outfile, 'w')
    
    # Write the header
    py.write("import numpy as np\n")
    # Write the function call
    py.write("def define_molec(num_spec, blength)\n")
    py.write('    mol_file = "tmp/%s.xyz"\n' % molname)
    py.write("    mol = open(mol_file, 'w')\n")
    py.write(r'    mol.write("%d\n")' % nchar[0])
    py.write('\n')
    py.write(r'    mol.write("%s\n")' % molname)
    py.write('\n')
    for i in range(nchar[0]):
        py.write(r'    mol.write("%s %.5f %.5f %.5f\n")' % (names[i],coords[i][1],coords[i][2],coords[i][3]))
        py.write('\n')
    py.write('    mol.close()\n')
    # Write Atoms Section
    py.write('    # Atomic Parameters\n')
    py.write('    atms = {\n')
    py.write('%s' % write_atm_line('name', names,nchar))
    py.write('%s' % write_atm_line('types', np.transpose(types)[1]-np.min(np.transpose(types)[1])+1,nchar))
    py.write('%s' % write_atm_line('q', np.transpose(charges)[1],nchar))
    py.write('%s' % write_atm_line('eps', names,nchar))
    py.write('%s' % write_atm_line('rmin', names,nchar))
    py.write('%s' % write_atm_line('sig', names,nchar))
    py.write('%s' % write_atm_line('mass', np.transpose(masses)[1],nchar))
    py.write('%s' % write_atm_line('group', zeros,nchar))
    py.write('%s' % write_atm_line('typgrp', zeros,nchar))
    py.write("        'cf':")
    py.write('["2 0 0"]\n')
    py.write('    }\n')
    
    py.write
    py.close()



if __name__ == '__main__':
    import sys
    import numpy as np

    # Read command line input
    if len(sys.argv) == 4:
        infile      = str(sys.argv[1])
        outfile     = str(sys.argv[2])
        molname     = str(sys.argv[3])
    elif len(sys.argv) == 2 and "-h" in sys.argv[1]:
        print("This is a program to create the python files to be read by the general build code.")
        print("This takes a connectivity file (as described https://lammps.sandia.gov/doc/molecule.html)")
        print("and a set of force field files (as read for lammps inputs)")
        print("in the form lmps.paircoeffs lmps.bondcoeffs etc")
        print("and this generates everything you need - reducing most of your work!")
        sys.exit("Good bye")
    else:
        print("Usage python molec_generator.py infile outfile molname")
        print("or python molec_generator.py -h for more help")
        sys.exit("Good bye")
    
    # Main function call
    write_molec_py(infile,outfile,molname)
    


