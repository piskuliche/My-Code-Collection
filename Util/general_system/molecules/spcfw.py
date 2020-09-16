import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/spcfw.xyz"
    mol = open(mol_file, 'w')
    mol.write("3\n")
    mol.write("spce\n")
    mol.write("O1 0.0000000 0.000000 0.0\n")
    mol.write("H1 1.0120000 0.000000 0.0 \n")
    mol.write("H2 -0.399318 0.929886 0.0\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["O1","H1","H2"],
        'atype':[1,2,2],
        'q':[-0.82,0.41,0.41],
        'eps':[0.1554253,0.0,0.0],
        'rmin':[],
        'sig':[3.165492,0.0,0.0],
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
       "O1-H1":[1, 528.5810000, 1.012, 0, 1],
        #O1-H2
       "O1-H2":[1, 528.5810000, 1.012, 0, 2]
    })
    # Angle Parameters
    angs = {
        'name':["H1-O1-H2"]
    }
    angs.update({
        #H1-O1-H2
       "H1-O1-H2":[1, 37.95000, 113.24000, 1, 0, 2]
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
