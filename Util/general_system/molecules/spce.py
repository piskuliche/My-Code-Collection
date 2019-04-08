import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/spce.xyz"
    mol = open(mol_file, 'w')
    mol.write("3\n")
    mol.write("spce\n")
    mol.write("O1 4.35567 2.73031 21.62297\n")
    mol.write("H1 4.77336 3.62526 21.77976\n")
    mol.write("H2 3.94618 2.39999 22.47339\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["O1","H1","H2"],
        'atype':[1,2,2],
        'q':[-0.8476,0.4238,0.4238],
        'eps':[0.15535,0.0,0.0],
        'rmin':[],
        'sig':[3.166,0.0,0.0],
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
       "O1-H1":[1, 553.00000, 1.00000, 0, 1],
        #O1-H2
       "O1-H2":[1, 553.00000, 1.00000, 0, 2]
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
