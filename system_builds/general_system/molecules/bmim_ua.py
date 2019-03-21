import numpy as np
def define_molec(num_spec, blength):
    # This is the bmim structure in my lab notebook on p. 46 with most h atoms removed.
    bmim_file = "tmp/bmim_ua.xyz"
    bmim = open(bmim_file, 'w')
    bmim.write("13\n")
    bmim.write("bmim_ua\n")
    bmim.write("N1    0.000   -0.647   -0.704\n")
    bmim.write("N3   0.000    1.400    0.000\n")
    bmim.write("C2    0.000    0.624   -1.068\n")
    bmim.write("H1   -0.000    0.975   -2.079\n")
    bmim.write("C5    0.000   -0.681    0.676\n")
    bmim.write("H3    0.000   -1.601    1.219\n")
    bmim.write("C4    0.000    0.582    1.110\n")
    bmim.write("H2    0.000    0.974    2.104\n")
    bmim.write("C6   -0.000    2.872    0.038\n")
    bmim.write("C7   -0.000   -1.861   -1.566\n")
    bmim.write("C8    0.000   -1.591   -3.080\n")
    bmim.write("C9    0.000   -2.870   -3.931\n")
    bmim.write("C10    0.000   -2.579   -5.434\n")
    bmim.close()
    # Atomic Parameters
    atms  = {
            #         0    1     2    3   4    5     6    7    8   9    10   11   12
            'name': ["N1","N3","C2","H1","C5","H3","C4","H2","C6","C7","C8","C9","C10"],
            'atype':[1, 1, 2, 3, 4, 5, 4, 5, 6, 7, 7, 7, 8 ],
            'q':[0.0454, -0.0067, -0.0794, 0.2512, -0.2072, 0.2544, -0.1384, 0.2295, 0.3235, 0.2519, -0.0034, 0.0406, 0.0386],
            'eps':[0.7113, 0.7113, 0.3598, 0.0628, 0.3598, 0.0628, 0.3598, 0.0628, 0.7894, 0.5962, 0.5580, 0.5580, 0.7728],
            'rmin':[],
            'sig':[3.250, 3.250, 3.400, 1.782, 3.400, 2.511, 3.400, 2.511, 3.813, 3.822, 3.947, 3.947, 3.902],
            'mass':[14.007, 14.007, 12.011, 1.008, 12.011, 1.008, 12.011, 1.008, 15.035, 14.027, 14.027, 14.027, 15.035],
            'group':[1, 1, 1, 2, 1, 2, 1, 2, 2, 3, 3, 3, 3],
            'typgrp':[0,0,0,1,0,1,0,1,1,1,1,1,1],
            'cf':["5 0 1 2 4 6", "4 3 5 7 8", "4 9 10 11 12"]
            }
    print("net charge is %s" % np.sum(atms['q']))
    # Bonded Parameters
    bnds = {
            'name':["C6-N3", "C7-N1", "C5-N1","C4-N3","C2-N1","C2-N3","C4-C5","C2-H1","C4-H2","C5-H3","C7-C8", "C8-C9", "C9-C10"]
            }
    bnds.update({
            # General Format Bond: Type, kb, ro, atm1, atm2
            # units: K: kJ/mol
            # units: r: Angstroms
            #CN3-NA
            "C6-N3":[1, 1409, 1.472,atms["name"].index(bnds['name'][0][:2]),atms["name"].index(bnds['name'][0][3:])],
            #CN2-NA
            "C7-N1":[1, 1409, 1.472,atms["name"].index(bnds['name'][1][:2]),atms["name"].index(bnds['name'][1][3:])],
            #CW-NA
            "C5-N1":[2, 1506, 1.378,atms["name"].index(bnds['name'][2][:2]),atms["name"].index(bnds['name'][2][3:])],
            "C4-N3":[2, 1506, 1.378,atms["name"].index(bnds['name'][3][:2]),atms["name"].index(bnds['name'][3][3:])],
            #CR-NA
            "C2-N1":[3, 1674, 1.325,atms["name"].index(bnds['name'][4][:2]),atms["name"].index(bnds['name'][4][3:])],
            "C2-N3":[3, 1674, 1.325,atms["name"].index(bnds['name'][5][:2]),atms["name"].index(bnds['name'][5][3:])],
            #CW-CW
            "C4-C5":[4, 1715, 1.343,atms["name"].index(bnds['name'][6][:2]),atms["name"].index(bnds['name'][6][3:])],
            #CR-H5
            "C2-H1":[5, 1590, 1.070,atms["name"].index(bnds['name'][7][:2]),atms["name"].index(bnds['name'][7][3:])],
            #CW-H4
            "C4-H2":[6, 1611, 1.070,atms["name"].index(bnds['name'][8][:2]),atms["name"].index(bnds['name'][8][3:])],
            "C5-H3":[6, 1611, 1.070,atms["name"].index(bnds['name'][9][:2]),atms["name"].index(bnds['name'][9][3:])],
            #CN2-CT2
            "C7-C8":[7, 1087, 1.526,atms["name"].index(bnds['name'][10][:2]),atms["name"].index(bnds['name'][10][3:])], 
            #CT2-CT2
            "C8-C9":[7, 1087, 1.526,atms["name"].index(bnds['name'][11][:2]),atms["name"].index(bnds['name'][11][3:])], 
            #CT2-CT3
            "C9-C10":[7, 1087, 1.526,atms["name"].index(bnds['name'][12][:2]),atms["name"].index(bnds['name'][12][3:])]
            })
    angs = {
            'name':["C8-C7-N1","C5-N1-C2","C4-N3-C2","N1-C2-H1","N3-C2-H1",
                    "H2-C4-C5","H3-C5-C4","N3-C4-H2","N1-C5-H3","C7-C8-C9","C8-C9-C10"]
            }
    angs.update({
            #CT2-CN2-Na
            "C8-C7-N1":[1,293,112.2,atms["name"].index(angs['name'][0][:2]),atms["name"].index(angs['name'][0][3:5]),atms["name"].index(angs['name'][0][6:])],
            #CW-Na-CR
            "C5-N1-C2":[2,502,108.0,atms["name"].index(angs['name'][1][:2]),atms["name"].index(angs['name'][1][3:5]),atms["name"].index(angs['name'][1][6:])],
            "C4-N3-C2":[2,502,108.0,atms["name"].index(angs['name'][2][:2]),atms["name"].index(angs['name'][2][3:5]),atms["name"].index(angs['name'][2][6:])],
            #NA-CR-H5
            "N1-C2-H1":[3,126,125.7,atms["name"].index(angs['name'][3][:2]),atms["name"].index(angs['name'][3][3:5]),atms["name"].index(angs['name'][3][6:])],
            "N3-C2-H1":[3,126,125.7,atms["name"].index(angs['name'][4][:2]),atms["name"].index(angs['name'][4][3:5]),atms["name"].index(angs['name'][4][6:])],
            #H4-CW-CW
            "H2-C4-C5":[4,126,130.7,atms["name"].index(angs['name'][5][:2]),atms["name"].index(angs['name'][5][3:5]),atms["name"].index(angs['name'][5][6:])],
            "H3-C5-C4":[4,126,130.7,atms["name"].index(angs['name'][6][:2]),atms["name"].index(angs['name'][6][3:5]),atms["name"].index(angs['name'][6][6:])],
            #NA-CW-H4
            "N3-C4-H2":[5,126,122.1,atms["name"].index(angs['name'][7][:2]),atms["name"].index(angs['name'][7][3:5]),atms["name"].index(angs['name'][7][6:])],
            "N1-C5-H3":[5,126,122.1,atms["name"].index(angs['name'][8][:2]),atms["name"].index(angs['name'][8][3:5]),atms["name"].index(angs['name'][8][6:])],
            #CN2-CT2-CT2
            "C7-C8-C9":[6,263,109.5,atms["name"].index(angs['name'][9][:2]),atms["name"].index(angs['name'][9][3:5]),atms["name"].index(angs['name'][9][6:])],
            #CT2-CT2-CT3
            "C8-C9-C10":[7,263,109.5,atms["name"].index(angs['name'][10][:2]),atms["name"].index(angs['name'][10][3:5]),atms["name"].index(angs['name'][10][6:])]
            })
    dihs = {
            'name':["C2-N3-C4-C5","N1-C5-C4-N3","N1-C2-N3-C4", "H1-C2-N1-C5","H1-C2-N3-C4", "H2-C4-C5-H3", "C4-C5-N1-C7", "C5-C4-N3-C6",
                    "H2-C4-N3-C2", "N1-C5-C4-H2", "N3-C4-C5-H3", "N1-C2-N3-C6", "N3-C2-N1-C7", "H1-C2-N3-C6", "H1-C2-N1-C7",
                    "H2-C4-N3-C6", "H3-C5-N1-C7", "C2-N1-C7-C8",
                    "C5-N1-C7-C8", "N1-C7-C8-C9", "C7-C8-C9-C10"]
            }
    dihs.update({
            #CR-NA-CW-CW
            "C2-N3-C4-C5":[1,50.21,2,180,atms["name"].index(dihs['name'][0][:2]),atms["name"].index(dihs['name'][0][3:5]),atms["name"].index(dihs['name'][0][6:8]),atms["name"].index(dihs['name'][0][9:])],
            #NA-CW-CW-NA
            "N1-C5-C4-N3":[1,50.21,2,180,atms["name"].index(dihs['name'][1][:2]),atms["name"].index(dihs['name'][1][3:5]),atms["name"].index(dihs['name'][1][6:8]),atms["name"].index(dihs['name'][1][9:])],
            #NA-CR-NA-CW
            "N1-C2-N3-C4":[1,50.21,2,180,atms["name"].index(dihs['name'][2][:2]),atms["name"].index(dihs['name'][2][3:5]),atms["name"].index(dihs['name'][2][6:8]),atms["name"].index(dihs['name'][2][9:])], 
            #H5-CR-NA-CW
            "H1-C2-N1-C5":[2,6.276,2,180,atms["name"].index(dihs['name'][3][:2]),atms["name"].index(dihs['name'][3][3:5]),atms["name"].index(dihs['name'][3][6:8]),atms["name"].index(dihs['name'][3][9:])],
            "H1-C2-N3-C4":[2,6.276,2,180,atms["name"].index(dihs['name'][4][:2]),atms["name"].index(dihs['name'][4][3:5]),atms["name"].index(dihs['name'][4][6:8]),atms["name"].index(dihs['name'][4][9:])], 
            #H4-CW-CW-H4
            "H2-C4-C5-H3":[2,6.276,2,180,atms["name"].index(dihs['name'][5][:2]),atms["name"].index(dihs['name'][5][3:5]),atms["name"].index(dihs['name'][5][6:8]),atms["name"].index(dihs['name'][5][9:])],
            #CW-CW-NA-CN2
            "C4-C5-N1-C7":[3,8.368,2,180,atms["name"].index(dihs['name'][6][:2]),atms["name"].index(dihs['name'][6][3:5]),atms["name"].index(dihs['name'][6][6:8]),atms["name"].index(dihs['name'][6][9:])], 
            #CW-CW-NA-CN3
            "C5-C4-N3-C6":[3,8.368,2,180,atms["name"].index(dihs['name'][7][:2]),atms["name"].index(dihs['name'][7][3:5]),atms["name"].index(dihs['name'][7][6:8]),atms["name"].index(dihs['name'][7][9:])],
            #H4-CW-NA-CR
            "H2-C4-N3-C2":[3,8.368,2,180,atms["name"].index(dihs['name'][8][:2]),atms["name"].index(dihs['name'][8][3:5]),atms["name"].index(dihs['name'][8][6:8]),atms["name"].index(dihs['name'][8][9:])], 
            #NA-CW-CW-H4
            "N1-C5-C4-H2":[2,6.276,2,180,atms["name"].index(dihs['name'][9][:2]),atms["name"].index(dihs['name'][9][3:5]),atms["name"].index(dihs['name'][9][6:8]),atms["name"].index(dihs['name'][9][9:])], 
            "N3-C4-C5-H3":[2,6.276,2,180,atms["name"].index(dihs['name'][10][:2]),atms["name"].index(dihs['name'][10][3:5]),atms["name"].index(dihs['name'][10][6:8]),atms["name"].index(dihs['name'][10][9:])],
            #NA-CR-NA-CN3
            "N1-C2-N3-C6":[3,8.368,2,180,atms["name"].index(dihs['name'][11][:2]),atms["name"].index(dihs['name'][11][3:5]),atms["name"].index(dihs['name'][11][6:8]),atms["name"].index(dihs['name'][11][9:])],
            #NA-CR-NA-CN2
            "N3-C2-N1-C7":[3,8.368,2,180,atms["name"].index(dihs['name'][12][:2]),atms["name"].index(dihs['name'][12][3:5]),atms["name"].index(dihs['name'][12][6:8]),atms["name"].index(dihs['name'][12][9:])],
            #H5-CR-NA-CN3
            "H1-C2-N3-C6":[2,6.276,2,180,atms["name"].index(dihs['name'][13][:2]),atms["name"].index(dihs['name'][13][3:5]),atms["name"].index(dihs['name'][13][6:8]),atms["name"].index(dihs['name'][13][9:])],
            #H5-CR-NA-CN2
            "H1-C2-N1-C7":[2,6.276,2,180,atms["name"].index(dihs['name'][14][:2]),atms["name"].index(dihs['name'][14][3:5]),atms["name"].index(dihs['name'][14][6:8]),atms["name"].index(dihs['name'][14][9:])],
            #H4-CW-NA-CN3
            "H2-C4-N3-C6":[2,6.276,2,180,atms["name"].index(dihs['name'][15][:2]),atms["name"].index(dihs['name'][15][3:5]),atms["name"].index(dihs['name'][15][6:8]),atms["name"].index(dihs['name'][15][9:])],
            #H4-CW-NA-CN2
            "H3-C5-N1-C7":[2,6.276,2,180,atms["name"].index(dihs['name'][16][:2]),atms["name"].index(dihs['name'][16][3:5]),atms["name"].index(dihs['name'][16][6:8]),atms["name"].index(dihs['name'][16][9:])],
            #CR-NA-CN2-CT2
            "C2-N1-C7-C8":[4,-0.987,1,0,atms["name"].index(dihs['name'][17][:2]),atms["name"].index(dihs['name'][17][3:5]),atms["name"].index(dihs['name'][17][6:8]),atms["name"].index(dihs['name'][17][9:])],
            #CW-NA-CN2-CT2
            "C5-N1-C7-C8":[5,-0.745,1,0,atms["name"].index(dihs['name'][18][:2]),atms["name"].index(dihs['name'][18][3:5]),atms["name"].index(dihs['name'][18][6:8]),atms["name"].index(dihs['name'][18][9:])],
            #NA-CN2-CT2-CT2
            "N1-C7-C8-C9":[6,4.180,3,0,atms["name"].index(dihs['name'][19][:2]),atms["name"].index(dihs['name'][19][3:5]),atms["name"].index(dihs['name'][19][6:8]),atms["name"].index(dihs['name'][19][9:])], 
            #CN2-CT2-CT2-CT3
            "C7-C8-C9-C10":[6,4.180,3,0,atms["name"].index(dihs['name'][20][:2]),atms["name"].index(dihs['name'][20][3:5]),atms["name"].index(dihs['name'][20][6:8]),atms["name"].index(dihs['name'][20][9:])]
            })
    imps = {
            'name':["H1-C2-N1-N3","C7-N1-C2-C5","C4-N3-C2-C6","H2-C4-C5-N3","H3-C5-N1-C4"]
            }
    imps.update({
            #NA-NA-CR-H5
            "H1-C2-N1-N3":[1,4.602,2,180,atms["name"].index(imps['name'][0][:2]),atms["name"].index(imps['name'][0][3:5]),atms["name"].index(imps['name'][0][6:8]),atms["name"].index(imps['name'][0][9:])],
            #CR-CW-NA-CN2           
            "C7-N1-C2-C5":[2,4.184,2,180,atms["name"].index(imps['name'][1][:2]),atms["name"].index(imps['name'][1][3:5]),atms["name"].index(imps['name'][1][6:8]),atms["name"].index(imps['name'][1][9:])],
            #CR-CW-NA-CN3
            "C4-N3-C2-C6":[2,4.184,2,180,atms["name"].index(imps['name'][2][:2]),atms["name"].index(imps['name'][2][3:5]),atms["name"].index(imps['name'][2][6:8]),atms["name"].index(imps['name'][2][9:])],
            #CW-NA-CW-H4
            "H2-C4-C5-N3":[1,4.602,2,180,atms["name"].index(imps['name'][3][:2]),atms["name"].index(imps['name'][3][3:5]),atms["name"].index(imps['name'][3][6:8]),atms["name"].index(imps['name'][3][9:])],
            "H3-C5-N1-C4":[1,4.602,2,180,atms["name"].index(imps['name'][4][:2]),atms["name"].index(imps['name'][4][3:5]),atms["name"].index(imps['name'][4][6:8]),atms["name"].index(imps['name'][4][9:])]
            })
    print("Note: The Liu bmim_ua forcefield has bond style of type:         K(r-req)^2")
    print("Note: The liu bmim_ua forcefield has angle style of type:        K(th-theq)^2")
    print("Note: The liu bmim_ua forcefield has dihedral style of type:     K(xi-xieq)^2")
    print("Note: The liu bmim_ua forcefield has improper style of type:     K(xi-xieq)^2")
    # Store the number of characteristics
    nchar = [13,13,11,21,5]
    ntyps = [11,8,7,6,2]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
