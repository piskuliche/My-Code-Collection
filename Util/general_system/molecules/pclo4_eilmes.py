import numpy as np
def define_molec(num_spec, blength):
    mol_file = "tmp/pclo4_eilmes.xyz"
    mol = open(mol_file, 'w')
    mol.write("5\n")
    mol.write("pclo4_eilmes\n")
    mol.write("Cl -0.00000 -0.00000 0.00000\n")
    mol.write("O1 1.54700 -0.00000 0.00000\n")
    mol.write("O2 -0.51570 1.45850 0.00000\n")
    mol.write("O3 -0.51570 -0.72930 1.26310\n")
    mol.write("O4 -0.51570 -0.72930 -1.26310\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["Cl","O1","O2","O3","O4"],
        'atype':[1,2,2,2,2],
        'q':[0.974754,-0.4936885,-0.4936885,-0.4936885,-0.4936885],
        'eps':[0.25,0.1561,0.1561,0.1561,0.1561],
        'rmin':[],
        'sig':[1.998,1.7959,1.7959,1.7959,1.7959],
        'mass':[35.453,15.999,15.999,15.999,15.999],
        'group':[0,0,0,0,0],
        'typgrp':[0,0,0,0,0],
        'cf':["2 0 0"]
    }
    # Bonding Parameters
    bnds = {
        'name':["Cl-O1","Cl-O2","Cl-O3","Cl-O4"]
    }
    bnds.update({
        # general format bond: type, kb, ro, atm1, tm2
        # Units kb:
        # Units r:
        #Cl-O1
       "Cl-O1":[1, 485.40000, 1.50600, 0, 1],
        #Cl-O2
       "Cl-O2":[1, 485.40000, 1.50600, 0, 2],
        #Cl-O3
       "Cl-O3":[1, 485.40000, 1.50600, 0, 3],
        #Cl-O4
       "Cl-O4":[1, 485.40000, 1.50600, 0, 4]
    })
    # Angle Parameters
    angs = {
        'name':["O1-Cl-O2","O1-Cl-O3","O1-Cl-O4","O2-Cl-O3","O2-Cl-O4","O3-Cl-O4"]
    }
    angs.update({
        #O1-Cl-O2
       "O1-Cl-O2":[1, 113.50000, 109.50000, 1, 0, 2],
        #O1-Cl-O3
       "O1-Cl-O3":[1, 113.50000, 109.50000, 1, 0, 3],
        #O1-Cl-O4
       "O1-Cl-O4":[1, 113.50000, 109.50000, 1, 0, 4],
        #O2-Cl-O3
       "O2-Cl-O3":[1, 113.50000, 109.50000, 2, 0, 3],
        #O2-Cl-O4
       "O2-Cl-O4":[1, 113.50000, 109.50000, 2, 0, 4],
        #O3-Cl-O4
       "O3-Cl-O4":[1, 113.50000, 109.50000, 3, 0, 4]
    })
    # Dihedral Parameters
    dihs = {
        'name':[]
    }
    # Improper Parameters
    imps = {
        'name':[]
    }
    nchar = [5, 4, 6, 0, 0]
    ntyps = [2, 1, 1, 0, 0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
