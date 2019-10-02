import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/acn_hirata.xyz"
    mol = open(mol_file, 'w')
    mol.write("3\n")
    mol.write("acn_maroncelli\n")
    mol.write("C1 -1.46000 0.00000 0.00000\n")
    mol.write("C2 0.00000 0.00000 0.00000\n")
    mol.write("N1 1.17000 0.00000 0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["C1","C2","N1"],
        'atype':[1,2,3],
        'q':[0.269,0.129,-0.398],
        'eps':[0.18004,0.2087,0.096975],
        'rmin':[],
        'sig':[3.8, 3.0, 3.4],
        'mass':[15.035,12.011,14.007],
        'group':[0,0,0],
        'typgrp':[0,0,0],
        'cf':["2 0 0"]
    }
    # Bonding Parameters
    bnds = {
        'name':["C1-C2","C2-N1"]
    }
    bnds.update({
        # general format bond: type, kb, ro, atm1, tm2
        # Units kb:
        # Units r:
        #O1-H1
       "C1-C2":[1, 190.00000, 1.46000, 0, 1],
        #H1-H2
       "C2-N1":[2, 640.00000, 1.17000, 1, 2]
    })
    # Angle Parameters
    angs = {
        'name':["C1-C2-N1"]
    }
    angs.update({
        #C1-C2-N1
       "C1-C2-N1":[1, 10.00000, 180.00000, 0, 1, 2]
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
