
import numpy as np
def define_molec(num_spec, blength):
    R32_file = "tmp/r32.xyz"
    R32 = open(R32_file, 'w')
    R32.write("5\n")
    R32.write("R32\n")
    R32.write("C     0.000    0.000   0.503\n")
    R32.write("F     0.000    1.110  -0.291\n")
    R32.write("F     0.000   -1.110  -0.291\n")
    R32.write("H    -0.908    0.000   1.107\n")
    R32.write("H     0.908    0.000   1.107\n")
    R32.close()
    atms = {
        #       0     1    2    3    4
        'name':["C1","F1","F2","H1","H2"],
        'atype':[1,2,2,3,3],
        'q':[0.43960, -0.26138, -0.26138, 0.04158, 0.04158 ],
        # K
        'eps':[54.6, 44.0, 44.0, 7.9, 7.9],
        'rmin':[],
        # A
        'sig':[3.15, 2.94, 2.94, 2.2293, 2.2293],
        'mass':[12.011, 18.997,18.997,1.008, 1.008],
        'group':[1,1,1,1,1],
        'typgrp':[1,1,1,1,1],
        'cf':["1 0", "4 1 2 3 4"]
        }
    print("net charge is %s" % np.sum(atms['q']))
    bnds = {
            'name':["C1-F1", "C1-F2", "C1-H1","C1-H2"]
            }
    bnds.update({
            #kj/mol/A^2
            "C1-F1":[1, 1544.61, 1.369,atms["name"].index(bnds['name'][0][:2]),atms["name"].index(bnds['name'][0][3:])],
            "C1-F2":[1, 1544.61, 1.369,atms["name"].index(bnds['name'][1][:2]),atms["name"].index(bnds['name'][1][3:])],
            "C1-H1":[2, 1472.89, 1.094,atms["name"].index(bnds['name'][2][:2]),atms["name"].index(bnds['name'][2][3:])], 
            "C1-H2":[2, 1472.89, 1.094,atms["name"].index(bnds['name'][3][:2]),atms["name"].index(bnds['name'][3][3:])]
            })
    angs = {
            'name':["H1-C1-H2", "F1-C1-F2", "F1-C1-H1", "F1-C1-H2", "F2-C1-H1", "F2-C1-H2"]
            }
    angs.update({
            #kj/mol/rad^2
            "H1-C1-H2":[1,146.54,113.6,atms["name"].index(angs['name'][0][:2]),atms["name"].index(angs['name'][0][3:5]),atms["name"].index(angs['name'][0][6:])],
            "F1-C1-F2":[2,367.61,108.7,atms["name"].index(angs['name'][1][:2]),atms["name"].index(angs['name'][1][3:5]),atms["name"].index(angs['name'][1][6:])],
            "F1-C1-H1":[3,249.92,108.6,atms["name"].index(angs['name'][2][:2]),atms["name"].index(angs['name'][2][3:5]),atms["name"].index(angs['name'][2][6:])],
            "F1-C1-H2":[3,249.92,108.6,atms["name"].index(angs['name'][3][:2]),atms["name"].index(angs['name'][3][3:5]),atms["name"].index(angs['name'][3][6:])],
            "F2-C1-H1":[3,249.92,108.6,atms["name"].index(angs['name'][4][:2]),atms["name"].index(angs['name'][4][3:5]),atms["name"].index(angs['name'][4][6:])],
            "F2-C1-H2":[3,249.92,108.6,atms["name"].index(angs['name'][5][:2]),atms["name"].index(angs['name'][5][3:5]),atms["name"].index(angs['name'][5][6:])]
            })
    dihs = {
            'name':[]
            }
    imps = {
            'name':[]
            }
    print("Note: The Raabe r32 forcefield has bond style of type: K(r-req)^2")
    print("Note: The Raabe r32 forcefield has angle style of type: K(th-theq)^2")
    # Write characteristics to file
    nchar = [5,4,6,0,0]
    ntyps = [3,2,3,0,0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
