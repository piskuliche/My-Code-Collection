import numpy as np
def define_molec(num_spec, blength)
    mol_file = "tmp/testmol.xyz"
    mol = open(mol_file, 'w')
    mol.write("5\n")
    mol.write("testmol\n")
    mol.write("C1 1.00000 1.00000 1.00000\n")
    mol.write("F1 2.00000 2.00000 2.00000\n")
    mol.write("F2 3.00000 3.00000 3.00000\n")
    mol.write("H1 4.00000 4.00000 4.00000\n")
    mol.write("H2 5.00000 5.00000 5.00000\n")
    mol.close()
    # Atomic Parameters
    atms = {
        'name':["C1","F1","F2","H1","H2"],
        'types':["1","2","2","3","3"],
        'q':["0.4396","-0.2614","-0.2614","0.0416","0.0416"],
        'eps':["C1","F1","F2","H1","H2"],
        'rmin':["C1","F1","F2","H1","H2"],
        'sig':["C1","F1","F2","H1","H2"],
        'mass':["12.011","18.997","18.997","1.008","1.008"],
        'group':["0.0","0.0","0.0","0.0","0.0"],
        'typgrp':["0.0","0.0","0.0","0.0","0.0"],
        'cf':["2 0 0"]
    }
