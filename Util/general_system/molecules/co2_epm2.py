import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/co2_epm2.xyz"
    mol = open(mol_file, 'w')
    mol.write("3\n")
    mol.write("co2_epm2\n")
    mol.write("O -1.14900 0.00000 0.00000\n")
    mol.write("C 0.00000 0.00000 0.00000\n")
    mol.write("O 1.14900 0.00000 0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["O","C","O"],
        'atype':[1,2,1],
        'q':[-0.3256,0.6512,-0.3256],
        'eps':[0.15998,0.05589,0.15998],
        'rmin':[],
        'sig':[3.033,2.757,3.033],
        'mass':[15.999,12.011,15.999],
        'group':[0,0,0],
        'typgrp':[0,0,0],
        'cf':["2 0 0"]
    }
    # Bonding Parameters
    bnds = {
        'name':["O-C","C-O"]
    }
    bnds.update({
        # general format bond: type, kb, ro, atm1, tm2
        # Units kb:
        # Units r:
        #O-C
       "O-C":[1, 575.00000, 1.14900, 0, 1],
        #C-O
       "C-O":[1, 575.00000, 1.14900, 1, 2]
    })
    # Angle Parameters
    angs = {
        'name':["O-C-O"]
    }
    angs.update({
        #O-C-O
       "O-C-O":[1, 29.00000, 180.00000, 0, 1, 2]
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
