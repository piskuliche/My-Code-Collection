import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/na_reduc ed.xyz"
    mol = open(mol_file, 'w')
    mol.write("1\n")
    mol.write("na_reduced\n")
    mol.write("Na 0.00000 0.00000 0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["Na"],
        'atype':[1],
        'q':[0.8],
        'eps':[0.016013],
        'rmin':[],
        'sig':[2.876],
        'mass':[22.990],
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
