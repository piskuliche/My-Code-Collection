import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/opc3.xyz"
    mol = open(mol_file, 'w')
    mol.write("3\n")
    mol.write("opc3\n")
    mol.write("O1 0.00000 0.00000  0.00000\n")
    mol.write("H1 -0.7992 0.56517  0.00000\n")
    mol.write("H2  0.7992 0.56517  0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["O1","H1","H2"],
        'atype':[1,2,2],
        'q':[-0.89517,0.447585,0.447585],
        'eps':[0.1634,0.0,0.0],
        'rmin':[],
        'sig':[3.17427,0.0,0.0],
        'mass':[15.999,1.008,1.008],
        'group':[0,0,0],
        'typgrp':[0,0,0],
        'cf':["2 0 0"]
    }
    # Bonding Parameters
    bnds = {
        'name':["O1-H1","O1-H2"]
    }
    bnds.update({
        # general format bond: type, kb, ro, atm1, tm2
        # Units kb:
        # Units r:
        #O1-H1
       "O1-H1":[1, 553.00000, 0.97888, 0, 1],
        #O1-H2
       "O1-H2":[1, 553.00000, 0.97888, 0, 2]
    })
    # Angle Parameters
    angs = {
        'name':["H1-O1-H2"]
    }
    angs.update({
        #H1-O1-H2
       "H1-O1-H2":[1, 55.00000, 109.47000, 1, 0, 2]
    })
    # Dihedral Parameters
    dihs = {
        'name':[]
    }
    # Improper Parameters
    imps = {
        'name':[]
    }
    nchar = [3, 2, 1, 0, 0]
    ntyps = [2, 1, 1, 0, 0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
