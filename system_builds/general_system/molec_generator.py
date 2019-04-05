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
    else:
        ntyps.append(0)
    if len(angles) > 0:
        ntyps.append(np.max(np.transpose(angles)[1]))
    else:
        ntyps.append(0)
    if len(dihedrals) > 0:
        ntyps.append(np.max(np.transpose(dihedrals)[1]))
    else:
        ntyps.append(0)
    if len(impropers) > 0:
        ntyps.append(np.max(np.transpose(impropers)[1]))
    else:
        ntyps.append(0)
    return ntyps

def write_atm_line(secname, array,nchar):
    line = ("        '%s':[" % secname)
    for i in range(nchar[0]):
        if i != 0:
            line += (',%s' % array[i])
        else:
            line += ('%s' % array[i])
    line += ("],\n")
    return line

def write_atm_line_str(secname, array,nchar):
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

def read_paircoeffs(types,nchar,molname):
    coeffs = []
    pcoeffs = []
    with open(molname+'.paircoeffs') as lmps:
        for line in lmps:
            coeffs.append([int(line.split()[1]), int(line.split()[2]), float(line.split()[3]), float(line.split()[4])])
    typmax = np.min(np.transpose(types)[1])
    for i in range(nchar[0]):
        pcoeffs.append(coeffs[types[i][1]- typmax])
    return pcoeffs

def read_bondcoeffs(bonds, nchar,molname):
    coeffs = []
    bcoeffs =[]
    with open(molname+'.bondcoeffs') as lmps:
        for line in lmps:
            coeffs.append([int(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
    typmax = int(np.min(np.transpose(coeffs)[0]))
    for i in range(nchar[1]):
        bcoeffs.append(coeffs[bonds[i][1]-typmax])
    return bcoeffs

def read_anglecoeffs(angles, nchar,molname):
    coeffs = []
    acoeffs =[]
    with open(molname+'.anglecoeffs') as lmps:
        for line in lmps:
            coeffs.append([int(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
    typmax = int(np.min(np.transpose(coeffs)[0]))
    for i in range(nchar[2]):
        acoeffs.append(coeffs[angles[i][1]-typmax])
    return acoeffs

def read_dihedralcoeffs(dihedrals, nchar,molname):
    coeffs = []
    dcoeffs =[]
    with open(molname+'.dihedralcoeffs') as lmps:
        for line in lmps:
            coeffs.append([int(line.split()[1]), float(line.split()[2]), float(line.split()[3]),float(line.split()[4])])
    typmax = int(np.min(np.transpose(coeffs)[0]))
    for i in range(nchar[3]):
        dcoeffs.append(coeffs[dihedrals[i][1]-typmax])
    return dcoeffs

def read_impropercoeffs(impropers, nchar,molname):
    coeffs = []
    icoeffs =[]
    with open(molname+'.impropercoeffs') as lmps:
        for line in lmps:
            coeffs.append([int(line.split()[1]), float(line.split()[2]), float(line.split()[3]),float(line.split()[4])])
    typmax = int(np.min(np.transpose(coeffs)[0]))
    for i in range(nchar[3]):
        icoeffs.append(coeffs[impropers[i][1]-typmax])
    return icoeffs

def write_molec_py(infile,outfile,molname):
    print("Begginning write to %s" % outfile)
    nchar, coords, types, charges, masses, bonds, angles, dihedrals, impropers = read_connect_file(infile)
    ntyps = calc_ntypes(types, bonds, angles, dihedrals, impropers)
    names = read_names(molname)
    # Coeffs
    pcoeffs = read_paircoeffs(types,nchar,molname)
    if nchar[1] != 0:
        bcoeffs = read_bondcoeffs(bonds, nchar,molname)
    if nchar[2] != 0:
        acoeffs = read_anglecoeffs(angles, nchar,molname)
    if nchar[3] != 0:
        dcoeffs = read_dihedralcoeffs(dihedrals, nchar,molname)
    if nchar[4] != 0:
        Icoeffs = read_impropercoeffs(impropers, nchar,molname)

    zeros = np.zeros(nchar[0])

    py = open(outfile, 'w')
    
    # Write the header
    py.write("import numpy as np\n")
    # Write the function call
    py.write("def define_molec(num_spec, blength):\n")
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
    py.write('%s' % write_atm_line_str('name', names,nchar))
    py.write('%s' % write_atm_line('atype', np.transpose(types)[1]-np.min(np.transpose(types)[1])+1,nchar))
    py.write('%s' % write_atm_line('q', np.transpose(charges)[1],nchar))
    py.write('%s' % write_atm_line('eps', np.transpose(pcoeffs)[2],nchar))
    py.write("        'rmin':")
    py.write('[],\n')
    py.write('%s' % write_atm_line('sig', np.transpose(pcoeffs)[3],nchar))
    py.write('%s' % write_atm_line('mass', np.transpose(masses)[1],nchar))
    py.write('%s' % write_atm_line('group', zeros.astype(int),nchar))
    py.write('%s' % write_atm_line('typgrp', zeros.astype(int),nchar))
    py.write("        'cf':")
    py.write('["2 0 0"]\n')
    py.write('    }\n')
    # Write Bonds Section
    py.write('    # Bonding Parameters\n')
    py.write('    bnds = {\n')
    py.write("        'name':[")
    count = 0
    for bond in bonds:
        if count == 0:
            py.write('"%s-%s"' % (names[bond[2]-1],names[bond[3]-1]))
        else:
            py.write(',"%s-%s"' % (names[bond[2]-1],names[bond[3]-1]))
        count += 1
    py.write(']\n')
    py.write('    }\n')
    py.write('    bnds.update({\n')
    py.write('        # general format bond: type, kb, ro, atm1, tm2\n')
    py.write('        # Units kb:\n')
    py.write('        # Units r:\n')
    count = 0
    if nchar[1] != 0:
        typmax = np.min(np.transpose(bcoeffs)[0])
    for bond in bonds:
       count += 1
       i = count-1
       py.write('        #%s-%s\n' % (names[bond[2]-1],names[bond[3]-1]))
       if count != len(bonds):
           py.write('       "%s-%s":[%d, %.5f, %.5f, %d, %d],\n'% (names[bond[2]-1],names[bond[3]-1],bcoeffs[i][0]-typmax+1, bcoeffs[i][1], bcoeffs[i][2],bond[2]-1, bond[3]-1))
       else:
           py.write('       "%s-%s":[%d, %.5f, %.5f, %d, %d]\n'% (names[bond[2]-1],names[bond[3]-1],bcoeffs[i][0]-typmax+1, bcoeffs[i][1], bcoeffs[i][2],bond[2]-1, bond[3]-1))
    py.write('    })\n')
    # Write Angles Section
    py.write('    # Angle Parameters\n')
    py.write('    angs = {\n')
    py.write("        'name':[")
    count = 0
    for angle in angles:
        if count == 0:
            py.write('"%s-%s-%s"' % (names[angle[2]-1],names[angle[3]-1],names[angle[4]-1]))
        else:
            py.write(',"%s-%s-%s"' % (names[angle[2]-1],names[angle[3]-1],names[angle[4]-1]))
        count += 1
    py.write(']\n')
    py.write('    }\n')
    if nchar[2] > 0:
        py.write('    angs.update({\n')
        count = 0
        typmax = np.min(np.transpose(acoeffs)[0])
        for angle in angles:
           count += 1
           i = count-1
           py.write('        #%s-%s-%s\n' % (names[angle[2]-1],names[angle[3]-1],names[angle[4]-1]))
           if count != len(angles):
               py.write('       "%s-%s-%s":[%d, %.5f, %.5f, %d, %d, %d],\n'% (names[angle[2]-1],names[angle[3]-1],names[angle[4]-1],acoeffs[i][0]-typmax+1, acoeffs[i][1], acoeffs[i][2],angle[2]-1, angle[3]-1,angle[4]-1))
           else:
               py.write('       "%s-%s-%s":[%d, %.5f, %.5f, %d, %d, %d]'% (names[angle[2]-1],names[angle[3]-1],names[angle[4]-1],acoeffs[i][0]-typmax+1, acoeffs[i][1], acoeffs[i][2],angle[2]-1, angle[3]-1,angle[4]-1))
        py.write('\n    })\n') 
   
    # Write Dihedral Section
    py.write('    # Dihedral Parameters\n')
    py.write('    dihs = {\n')
    py.write("        'name':[")
    count = 0
    for dihedral in dihedrals:
        if count == 0:
            py.write('"%s-%s-%s-%s"' % (names[dihedral[2]-1],names[dihedral[3]-1],names[dihedral[4]-1],names[dihedral[5]-1]))
        else:
            py.write(',"%s-%s-%s-%s"' % (names[dihedral[2]-1],names[dihedral[3]-1],names[dihedral[4]-1],names[dihedral[5]-1]))
        count += 1
    py.write(']\n')
    py.write('    }\n')
    if nchar[3] > 0:
        py.write('    dihs.update({\n')
        count = 0
        typmax = np.min(np.transpose(dcoeffs)[0])
        for dihedral in dihedrals:
           count += 1
           i = count-1
           py.write('        #%s-%s-%s-%s\n' % (names[dihedral[2]-1],names[dihedral[3]-1],names[dihedral[4]-1],names[dihedral[5]-1]))
           if count != len(angles):
               py.write('       "%s-%s-%s-%s":[%d, %.5f, %d,%.5f, %d, %d, %d, %d],\n'% (names[dihedral[2]-1],names[dihedral[3]-1],names[dihedral[4]-1],names[dihedral[5]-1],dcoeffs[i][0]-typmax+1, dcoeffs[i][1], dcoeffs[i][2],dcoeffs[i][3],dihedral[2]-1, dihedral[3]-1,dihedral[4]-1,dihedral[5]-1))
           else:
               py.write('       "%s-%s-%s-%s":[%d, %.5f, %d,%.5f, %d, %d, %d, %d]'% (names[dihedral[2]-1],names[dihedral[3]-1],names[dihedral[4]-1],names[dihedral[5]-1],dcoeffs[i][0]-typmax+1, dcoeffs[i][1], dcoeffs[i][2],dcoeffs[i][3],dihedral[2]-1, dihedral[3]-1,dihedral[4]-1,dihedral[5]-1))
        py.write('\n    })\n')
    # Write Improper Section
    py.write('    # Improper Parameters\n')
    py.write('    imps = {\n')
    py.write("        'name':[")
    count = 0
    for improper in impropers:
        if count == 0:
            py.write('"%s-%s-%s-%s"' % (names[improper[2]-1],names[improper[3]-1],names[improper[4]-1],names[improper[5]-1]))
        else:
            py.write(',"%s-%s-%s-%s"' % (names[improper[2]-1],names[improper[3]-1],names[improper[4]-1],names[improper[5]-1]))
        count += 1
    py.write(']\n')
    py.write('    }\n')
    if nchar[3] > 0:
        py.write('    imps.update({\n')
        count = 0
        typmax = np.min(np.transpose(icoeffs)[0])
        for improper in impropers:
           count += 1
           i = count-1
           py.write('        #%s-%s-%s-%s\n' % (names[improper[2]-1],names[improper[3]-1],names[improper[4]-1],names[improper[5]-1]))
           if count != len(angles):
               py.write('       "%s-%s-%s-%s":[%d, %.5f, %d,%.5f, %d, %d, %d, %d],\n'% (names[improper[2]-1],names[improper[3]-1],names[improper[4]-1],names[improper[5]-1],icoeffs[i][0]-typmax+1, icoeffs[i][1], icoeffs[i][2],icoeffs[i][3],improper[2]-1, improper[3]-1,improper[4]-1,improper[5]-1))
           else:
               py.write('       "%s-%s-%s-%s":[%d, %.5f, %d,%.5f, %d, %d, %d, %d]'% (names[improper[2]-1],names[improper[3]-1],names[improper[4]-1],names[improper[5]-1],icoeffs[i][0]-typmax+1, icoeffs[i][1], icoeffs[i][2],icoeffs[i][3],improper[2]-1, improper[3]-1,improper[4]-1,improper[5]-1))
        py.write('\n    })\n')
    
    py.write('    nchar = [%d, %d, %d, %d, %d]\n' % (nchar[0],nchar[1],nchar[2],nchar[3],nchar[4]))
    py.write('    ntyps = [%d, %d, %d, %d, %d]\n' % (ntyps[0],ntyps[1],ntyps[2],ntyps[3],ntyps[4]))
    py.write('    return nchar, ntyps, atms, bnds, angs, dihs, imps\n')
    py.close()
    print("%s write complete." % outfile)
    return



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
    print("Program Complete. Good bye.")
    


