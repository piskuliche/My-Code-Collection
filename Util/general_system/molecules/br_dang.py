import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/br_dang.xyz"
    mol = open(mol_file, 'w')
    mol.write("1\n")
    mol.write("br_dang\n")
    mol.write("Br 0.00000 0.00000 0.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["Br"],
        'atype':[1],
        'q':[-1.0],
        'eps':[0.12467],
        'rmin':[],
        'sig':[3.896],
        'mass':[79.904],
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
