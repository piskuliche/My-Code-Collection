import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/li_eilmes.xyz"
    mol = open(mol_file, 'w')
    mol.write("1\n")
    mol.write("li_eilmes\n")
    mol.write("Li 0.00000 0.00000 0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["Li"],
        'atype':[1],
        'q':[1.0],
        'eps':[0.018279],
        'rmin':[],
        'sig':[1.1371],
        'mass':[6.941],
        'group':[0],
        'typgrp':[0],
        'cf':["2 0 0"]
    }
    # Bonding Parameters
    bnds = {
        'name':[]
    }
    bnds.update({
        # general format bond: type, kb, ro, atm1, tm2
        # Units kb:
        # Units r:
    })
    # Angle Parameters
    angs = {
        'name':[]
    }
    # Dihedral Parameters
    dihs = {
        'name':[]
    }
    # Improper Parameters
    imps = {
        'name':[]
    }
    nchar = [1, 0, 0, 0, 0]
    ntyps = [1, 0, 0, 0, 0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
