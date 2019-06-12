import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/acn_maroncelli.xyz"
    mol = open(mol_file, 'w')
    mol.write("3\n")
    mol.write("acn_maroncelli\n")
    mol.write("O1 -1.46000 0.00000 0.00000\n")
    mol.write("H1 0.00000 0.00000 0.00000\n")
    mol.write("H2 1.17000 0.00000 0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["O1","H1","H2"],
        'atype':[1,2,3],
        'q':[0.269,0.129,-0.398],
        'eps':[1.588,0.4157,0.4157],
        'rmin':[],
        'sig':[3.6,3.4,3.3],
        'mass':[15.035,12.011,14.007],
        'group':[0,0,0],
        'typgrp':[0,0,0],
        'cf':["2 0 0"]
    }
    # Bonding Parameters
    bnds = {
        'name':["O1-H1","H1-H2"]
    }
    bnds.update({
        # general format bond: type, kb, ro, atm1, tm2
        # Units kb:
        # Units r:
        #O1-H1
       "O1-H1":[1, 190.00000, 1.46000, 0, 1],
        #H1-H2
       "H1-H2":[2, 640.00000, 1.17000, 1, 2]
    })
    # Angle Parameters
    angs = {
        'name':["O1-H1-H2"]
    }
    angs.update({
        #O1-H1-H2
       "O1-H1-H2":[1, 10.00000, 180.00000, 0, 1, 2]
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
    ntyps = [3, 2, 1, 0, 0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
