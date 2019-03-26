import numpy as np
def define_molec(num_spec, blength):
    # This is the PF6 structure in my lab notebook on p.44
    HFP_file = "tmp/pf6_ua.xyz"
    HFP = open(HFP_file, 'w')
    HFP.write("7\n")
    HFP.write("PF6\n")
    HFP.write("P     0.000   0.000   0.000\n")
    HFP.write("F     0.000   0.000   1.646\n")
    HFP.write("F     0.000   1.646   0.000\n")
    HFP.write("F     0.000   0.000  -1.646\n")
    HFP.write("F     0.000  -1.646   0.000\n")
    HFP.write("F    -1.646   0.000   0.000\n")
    HFP.write("F     1.646   0.000   0.000\n")
    HFP.close()
    atms = {
            #        0    1    2   3      4    5   6
            'name':["P1","F1","F2","F3","F4","F5","F6"],
            'atype':[1,2,2,2,2,2,2],
            'q':[0.7562, -0.2927, -0.2927, -0.2927, -0.2927, -0.2927, -0.2927],
            # Units KJ/mol
            'eps':[0.8368, 0.2552, 0.2552,0.2552,0.2552,0.2552,0.2552],
            'rmin':[],
            'sig':[3.742,3.118,3.118,3.118,3.118,3.118,3.118],
            'mass':[30.974,  18.997,18.997,18.997,18.997,18.997,18.997],
            'group':[1,1,1,1,1,1,1],
            'typgrp':[1,1,1,1,1,1,1],
            'cf':["1 0", "6 1 2 3 4 5 6"]
            }
    print("net charge is %s" % np.sum(atms['q']))
    bnds = {
            'name':["P1-F1", "P1-F2", "P1-F3","P1-F4","P1-F5","P1-F6"]
            }
    bnds.update({
            # units kj/mol/Ang^2
            "P1-F1":[1, 795.0, 1.600,atms["name"].index(bnds['name'][0][:2]),atms["name"].index(bnds['name'][0][3:])],
            "P1-F2":[1, 795.0, 1.600,atms["name"].index(bnds['name'][1][:2]),atms["name"].index(bnds['name'][1][3:])],
            "P1-F3":[1, 795.0, 1.600,atms["name"].index(bnds['name'][2][:2]),atms["name"].index(bnds['name'][2][3:])],
            "P1-F4":[1, 795.0, 1.600,atms["name"].index(bnds['name'][3][:2]),atms["name"].index(bnds['name'][3][3:])],
            "P1-F5":[1, 795.0, 1.600,atms["name"].index(bnds['name'][4][:2]),atms["name"].index(bnds['name'][4][3:])],
            "P1-F6":[1, 795.0, 1.600,atms["name"].index(bnds['name'][5][:2]),atms["name"].index(bnds['name'][5][3:])]
            })
    angs = {
            'name':["F1-P1-F2", "F1-P1-F6", "F1-P1-F4", "F1-P1-F5", "F2-P1-F5", "F5-P1-F4", "F4-P1-F6","F6-P1-F2", "F3-P1-F2", "F3-P1-F4", "F3-P1-F5","F3-P1-F6"]
            }
    angs.update({
            # Units kj/mol/rad^2
            "F1-P1-F2":[1,334.7,90.0,atms["name"].index(angs['name'][0][:2]),atms["name"].index(angs['name'][0][3:5]),atms["name"].index(angs['name'][0][6:])],
            "F1-P1-F6":[1,334.7,90.0,atms["name"].index(angs['name'][1][:2]),atms["name"].index(angs['name'][1][3:5]),atms["name"].index(angs['name'][1][6:])],
            "F1-P1-F4":[1,334.7,90.0,atms["name"].index(angs['name'][2][:2]),atms["name"].index(angs['name'][2][3:5]),atms["name"].index(angs['name'][2][6:])],
            "F1-P1-F5":[1,334.7,90.0,atms["name"].index(angs['name'][3][:2]),atms["name"].index(angs['name'][3][3:5]),atms["name"].index(angs['name'][3][6:])],
            "F2-P1-F5":[1,334.7,90.0,atms["name"].index(angs['name'][4][:2]),atms["name"].index(angs['name'][4][3:5]),atms["name"].index(angs['name'][4][6:])],
            "F5-P1-F4":[1,334.7,90.0,atms["name"].index(angs['name'][5][:2]),atms["name"].index(angs['name'][5][3:5]),atms["name"].index(angs['name'][5][6:])],
            "F4-P1-F6":[1,334.7,90.0,atms["name"].index(angs['name'][6][:2]),atms["name"].index(angs['name'][6][3:5]),atms["name"].index(angs['name'][6][6:])],
            "F6-P1-F2":[1,334.7,90.0,atms["name"].index(angs['name'][7][:2]),atms["name"].index(angs['name'][7][3:5]),atms["name"].index(angs['name'][7][6:])],
            "F3-P1-F2":[1,334.7,90.0,atms["name"].index(angs['name'][8][:2]),atms["name"].index(angs['name'][8][3:5]),atms["name"].index(angs['name'][8][6:])],
            "F3-P1-F4":[1,334.7,90.0,atms["name"].index(angs['name'][9][:2]),atms["name"].index(angs['name'][9][3:5]),atms["name"].index(angs['name'][9][6:])],
            "F3-P1-F5":[1,334.7,90.0,atms["name"].index(angs['name'][10][:2]),atms["name"].index(angs['name'][10][3:5]),atms["name"].index(angs['name'][10][6:])],
            "F3-P1-F6":[1,334.7,90.0,atms["name"].index(angs['name'][11][:2]),atms["name"].index(angs['name'][11][3:5]),atms["name"].index(angs['name'][11][6:])],
            })
    dihs = {
            'name':[]
            }
    imps = {
            'name':[]
            }
    print("Note: The liu pf6_ua forcefield has bond style of type: K(r-req)^2")
    print("Note: The liu pf6_ua forcefield has angle style of type: K(th-theq)^2")
    # Stores the number of characteristics
    nchar = [7,6,12,0,0]
    ntyps =[2,1,1,0,0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps

