def define_molec(num_spec, blength):
    # This is the bmim structure in my lab notebook on p. 44
    bmim_file = "tmp/bmim.xyz"
    bmim = open(bmim_file, 'w')
    bmim.write("25\n")
    bmim.write("bmim\n")
    bmim.write("N1    0.000   -0.647   -0.704\n")
    bmim.write("N3   0.000    1.400    0.000\n")
    bmim.write("C2    0.000    0.624   -1.068\n")
    bmim.write("H1   -0.000    0.975   -2.079\n")
    bmim.write("C5    0.000   -0.681    0.676\n")
    bmim.write("H3    0.000   -1.601    1.219\n")
    bmim.write("C4    0.000    0.582    1.110\n")
    bmim.write("H2    0.000    0.974    2.104\n")
    bmim.write("C6   -0.000    2.872    0.038\n")
    bmim.write("H6   -0.879    3.217    0.567\n")
    bmim.write("H4    0.000    3.293   -0.957\n")
    bmim.write("H5    0.879    3.217    0.567\n")
    bmim.write("C7   -0.000   -1.861   -1.566\n")
    bmim.write("H8   -0.896   -2.421   -1.332\n")
    bmim.write("H7    0.896   -2.421   -1.332\n")
    bmim.write("C8    0.000   -1.591   -3.080\n")
    bmim.write("H9    0.876   -0.999   -3.339\n")
    bmim.write("H10   -0.876   -0.999   -3.339\n")
    bmim.write("C9    0.000   -2.870   -3.931\n")
    bmim.write("H11   -0.871   -3.477   -3.688\n")
    bmim.write("H12    0.871   -3.477   -3.688\n")
    bmim.write("C10    0.000   -2.579   -5.434\n")
    bmim.write("H13    0.877   -2.017   -5.736\n")
    bmim.write("H14   -0.8777  -2.017   -5.736\n")
    bmim.write("H15    0.000   -3.501   -6.007\n")
    bmim.close()
    # Atomic Parameters
    atms  = {
            'name': ["N1","N3","C2","H1","C5","H3","C4","H2","C6","H6","H4","H5","C7", "H8", "H7", "C8","H9","H10","C9","H11","H12","C10","H13","H14","H15"],
            'atype':[1, 1, 2, 3, 2, 4, 2, 4, 5, 6, 6, 6, 5, 7, 7, 8, 9, 9, 8, 9, 9, 10, 11, 11, 11],
            'q':[0.111, 0.133, 0.056, 0.177, -0.217, 0.207, -0.141, 0.181, -0.157, 0.142, 0.152, 0.073, 0.095, 0.045, 0.055, -0.122, 0.055, 0.001, 0.256, -0.029, -0.099, -0.209, 0.051, 0.040, 0.075],
            'kjeps':[0.837, 0.837, 0.209, 0.192,0.209, 0.033,0.209, 0.033, 0.084, 0.092, 0.092, 0.092, 0.084, 0.092, 0.092, 0.234, 0.117, 0.117, 0.234,0.117, 0.117, 0.326, 0.100, 0.100, 0.100],
            'rmin':[3.700, 3.700, 3.600, 1.800,3.600,2.936, 3.600, 2.936, 4.550,2.640, 2.640, 2.640,4.550, 2.680, 2.680, 4.020, 2.680, 2.680,4.020,2.680,2.680, 4.080,2.680,2.680, 2.680],
            'mass':[14.007, 14.007, 12.011, 1.008, 12.011, 1.008, 12.011, 1.008,12.011, 1.008,1.008, 1.008, 12.011, 1.008, 1.008, 12.011, 1.008, 1.008, 12.011, 1.008, 1.008, 12.011, 1.008, 1.008, 1.008]
            }
    # Bonded Parameters
    bnds = {
            'name':["C6-N3", "C7-N1", "C5-N1","C4-N3","C2-N1","C2-N3","C4-C5","C2-H1","C4-H2","C5-H3","C6-H4","C6-H5","C6-H6","C7-H7",
                    "C7-H8","C8-H9","C8-H10","C9-H11","C9-H12","C10-H13","C10-H14","C10-H15","C7-C8", "C8-C9", "C9-C10"]
            }
    bnds.update({
            # General Format Bond: Type, kb, ro, atm1, atm2
            "C6-N3":[1, 220.0, 1.470,atms["name"].index(bnds['name'][0][:2]),atms["name"].index(bnds['name'][0][3:])],
            "C7-N1":[2, 220.0, 1.483,atms["name"].index(bnds['name'][1][:2]),atms["name"].index(bnds['name'][1][3:])],
            "C5-N1":[3, 400.0, 1.382,atms["name"].index(bnds['name'][2][:2]),atms["name"].index(bnds['name'][2][3:])],
            "C4-N3":[3, 400.0, 1.382,atms["name"].index(bnds['name'][3][:2]),atms["name"].index(bnds['name'][3][3:])],
            "C2-N1":[4, 400.0, 1.337,atms["name"].index(bnds['name'][4][:2]),atms["name"].index(bnds['name'][4][3:])],
            "C2-N3":[5, 400.0, 1.337,atms["name"].index(bnds['name'][5][:2]),atms["name"].index(bnds['name'][5][3:])],
            "C4-C5":[6, 410.0, 1.361,atms["name"].index(bnds['name'][6][:2]),atms["name"].index(bnds['name'][6][3:])],
            "C2-H1":[7, 0.0, 1.078,atms["name"].index(bnds['name'][7][:2]),atms["name"].index(bnds['name'][7][3:])],
            "C4-H2":[7, 0.0, 1.078,atms["name"].index(bnds['name'][8][:2]),atms["name"].index(bnds['name'][8][3:])],
            "C5-H3":[7, 0.0, 1.078,atms["name"].index(bnds['name'][9][:2]),atms["name"].index(bnds['name'][9][3:])],
            "C6-H4":[8, 0.0, 1.089,atms["name"].index(bnds['name'][10][:2]),atms["name"].index(bnds['name'][10][3:])],
            "C6-H5":[8, 0.0, 1.089,atms["name"].index(bnds['name'][11][:2]),atms["name"].index(bnds['name'][11][3:])],
            "C6-H6":[8, 0.0, 1.089,atms["name"].index(bnds['name'][12][:2]),atms["name"].index(bnds['name'][12][3:])],
            "C7-H7":[9, 0.0, 1.091,atms["name"].index(bnds['name'][13][:2]),atms["name"].index(bnds['name'][13][3:])],
            "C7-H8":[9, 0.0, 1.091,atms["name"].index(bnds['name'][14][:2]),atms["name"].index(bnds['name'][14][3:])],
            "C8-H9":[10, 0.0, 1.096,atms["name"].index(bnds['name'][15][:2]),atms["name"].index(bnds['name'][15][3:])],
            "C8-H10":[10, 0.0, 1.096,atms["name"].index(bnds['name'][16][:2]),atms["name"].index(bnds['name'][16][3:])],
            "C9-H11":[10, 0.0, 1.096,atms["name"].index(bnds['name'][17][:2]),atms["name"].index(bnds['name'][17][3:])],
            "C9-H12":[10, 0.0, 1.096,atms["name"].index(bnds['name'][18][:2]),atms["name"].index(bnds['name'][18][3:])],
            "C10-H13":[11, 0.0, 1.093,atms["name"].index(bnds['name'][19][:3]),atms["name"].index(bnds['name'][19][4:])],
            "C10-H14":[11, 0.0, 1.093,atms["name"].index(bnds['name'][20][:3]),atms["name"].index(bnds['name'][20][4:])],
            "C10-H15":[11, 0.0, 1.093,atms["name"].index(bnds['name'][21][:3]),atms["name"].index(bnds['name'][21][4:])],
            "C7-C8":[12, 200.0, 1.530,atms["name"].index(bnds['name'][22][:2]),atms["name"].index(bnds['name'][22][3:])],
            "C8-C9":[13, 222.5, 1.534,atms["name"].index(bnds['name'][23][:2]),atms["name"].index(bnds['name'][23][3:])],
            "C9-C10":[14,222.5, 1.646,atms["name"].index(bnds['name'][24][:2]),atms["name"].index(bnds['name'][24][3:])]
            })
    angs = {
            'name':["C8-C7-N1","C5-N1-C2","C4-N3-C2","H4-C6-N3","H5-C6-N3","H6-C6-N3","N1-C7-H8","N1-C7-H7","N1-C2-H1","N3-C2-H1",
                    "H2-C4-C5","H3-C5-C4","N3-C4-H2","H4-C6-H5","H4-C6-H6","H5-C6-H6","H7-C7-H8","H7-C7-C8","H8-C7-C8","H9-C8-C7",
                    "H10-C8-C7", "C7-C8-C9","C8-C9-C10","H9-C8-H10","H11-C9-H12","C8-C9-H11","C8-C9-H12","C9-C10-H13","C9-C10-H14",
                    "C9-C10-H15","H13-C10-H14","H13-C10-H15","H14-C10-H15","N1-C2-H1","N3-C2-H1"]
            }
    angs.update({
            "C8-C7-N1":[1,140.0,112.6,atms["name"].index(angs['name'][0][:2]),atms["name"].index(angs['name'][0][3:5]),atms["name"].index(angs['name'][0][6:])],
            "C5-N1-C2":[2,130.0,125.9,atms["name"].index(angs['name'][1][:2]),atms["name"].index(angs['name'][1][3:5]),atms["name"].index(angs['name'][1][6:])],
            "C4-N3-C2":[2,130.0,125.9,atms["name"].index(angs['name'][2][:2]),atms["name"].index(angs['name'][2][3:5]),atms["name"].index(angs['name'][2][6:])],
            "H4-C6-N3":[3,30.0,109.6,atms["name"].index(angs['name'][3][:2]),atms["name"].index(angs['name'][3][3:5]),atms["name"].index(angs['name'][3][6:])],
            "H5-C6-N3":[3,30.0,109.6,atms["name"].index(angs['name'][4][:2]),atms["name"].index(angs['name'][4][3:5]),atms["name"].index(angs['name'][4][6:])],
            "H6-C6-N3":[3,30.0,109.6,atms["name"].index(angs['name'][5][:2]),atms["name"].index(angs['name'][5][3:5]),atms["name"].index(angs['name'][5][6:])],
            "N1-C7-H8":[4,30.0,106.8,atms["name"].index(angs['name'][6][:2]),atms["name"].index(angs['name'][6][3:5]),atms["name"].index(angs['name'][6][6:])],
            "N1-C7-H7":[4,30.0,106.8,atms["name"].index(angs['name'][7][:2]),atms["name"].index(angs['name'][7][3:5]),atms["name"].index(angs['name'][7][6:])],
            "N1-C2-H1":[5,130.0,107.2,atms["name"].index(angs['name'][8][:2]),atms["name"].index(angs['name'][8][3:5]),atms["name"].index(angs['name'][8][6:])],
            "N3-C2-H1":[6,130.0,109.1,atms["name"].index(angs['name'][9][:2]),atms["name"].index(angs['name'][9][3:5]),atms["name"].index(angs['name'][9][6:])],
            "H2-C4-C5":[7,25.0,130.8,atms["name"].index(angs['name'][10][:2]),atms["name"].index(angs['name'][10][3:5]),atms["name"].index(angs['name'][10][6:])],
            "H3-C5-C4":[7,25.0,130.8,atms["name"].index(angs['name'][11][:2]),atms["name"].index(angs['name'][11][3:5]),atms["name"].index(angs['name'][11][6:])],
            "N3-C4-H2":[8,25.0,112.6,atms["name"].index(angs['name'][12][:2]),atms["name"].index(angs['name'][12][3:5]),atms["name"].index(angs['name'][12][6:])],
            "H4-C6-H5":[9,35.5,109.3,atms["name"].index(angs['name'][13][:2]),atms["name"].index(angs['name'][13][3:5]),atms["name"].index(angs['name'][13][6:])],
            "H4-C6-H6":[9,35.5,109.3,atms["name"].index(angs['name'][14][:2]),atms["name"].index(angs['name'][14][3:5]),atms["name"].index(angs['name'][14][6:])],
            "H5-C6-H6":[9,35.5,109.3,atms["name"].index(angs['name'][15][:2]),atms["name"].index(angs['name'][15][3:5]),atms["name"].index(angs['name'][15][6:])],
            "H7-C7-H8":[10,35.5,107.2,atms["name"].index(angs['name'][16][:2]),atms["name"].index(angs['name'][16][3:5]),atms["name"].index(angs['name'][16][6:])],
            "H7-C7-C8":[11,35.5,111.5,atms["name"].index(angs['name'][17][:2]),atms["name"].index(angs['name'][17][3:5]),atms["name"].index(angs['name'][17][6:])],
            "H8-C7-C8":[11,35.5,111.5,atms["name"].index(angs['name'][18][:2]),atms["name"].index(angs['name'][18][3:5]),atms["name"].index(angs['name'][18][6:])],
            "H9-C8-C7":[12,33.4,109.5,atms["name"].index(angs['name'][19][:2]),atms["name"].index(angs['name'][19][3:5]),atms["name"].index(angs['name'][19][6:])],
            "H10-C8-C7":[12,33.4,109.5,atms["name"].index(angs['name'][20][:3]),atms["name"].index(angs['name'][20][4:6]),atms["name"].index(angs['name'][20][7:])],
            "C7-C8-C9":[13,58.4,111.6,atms["name"].index(angs['name'][21][:2]),atms["name"].index(angs['name'][21][3:5]),atms["name"].index(angs['name'][21][6:])],
            "C8-C9-C10":[13,58.4,111.6,atms["name"].index(angs['name'][22][:2]),atms["name"].index(angs['name'][22][3:5]),atms["name"].index(angs['name'][22][6:])],
            "H9-C8-H10":[14,34.5,106.4,atms["name"].index(angs['name'][23][:2]),atms["name"].index(angs['name'][23][3:5]),atms["name"].index(angs['name'][23][6:])],
            "H11-C9-H12":[14,34.5,106.4,atms["name"].index(angs['name'][24][:3]),atms["name"].index(angs['name'][24][4:6]),atms["name"].index(angs['name'][24][7:])],
            "C8-C9-H11":[15,34.6,109.7,atms["name"].index(angs['name'][25][:2]),atms["name"].index(angs['name'][25][3:5]),atms["name"].index(angs['name'][25][6:])],
            "C8-C9-H12":[15,34.6,109.7,atms["name"].index(angs['name'][26][:2]),atms["name"].index(angs['name'][26][3:5]),atms["name"].index(angs['name'][26][6:])],
            "C9-C10-H13":[15,34.6,109.7,atms["name"].index(angs['name'][27][:2]),atms["name"].index(angs['name'][27][3:6]),atms["name"].index(angs['name'][27][7:])],
            "C9-C10-H14":[15,34.6,109.7,atms["name"].index(angs['name'][28][:2]),atms["name"].index(angs['name'][28][3:6]),atms["name"].index(angs['name'][28][7:])],
            "C9-C10-H15":[15,34.6,109.7,atms["name"].index(angs['name'][29][:2]),atms["name"].index(angs['name'][29][3:6]),atms["name"].index(angs['name'][29][7:])],
            "H13-C10-H14":[16,35.5,107.6,atms["name"].index(angs['name'][30][:3]),atms["name"].index(angs['name'][30][4:7]),atms["name"].index(angs['name'][30][8:])],
            "H13-C10-H15":[16,35.5,107.6,atms["name"].index(angs['name'][31][:3]),atms["name"].index(angs['name'][31][4:7]),atms["name"].index(angs['name'][31][8:])],
            "H14-C10-H15":[16,35.5,107.6,atms["name"].index(angs['name'][32][:3]),atms["name"].index(angs['name'][32][4:7]),atms["name"].index(angs['name'][32][8:])],
            "N1-C2-H1":[17,25.0,125.5,atms["name"].index(angs['name'][33][:2]),atms["name"].index(angs['name'][33][3:5]),atms["name"].index(angs['name'][33][6:])],
            "N3-C2-H1":[17,25.0,125.5,atms["name"].index(angs['name'][34][:2]),atms["name"].index(angs['name'][34][3:5]),atms["name"].index(angs['name'][34][6:])]
            })
    dihs = {
            'name':["C2-N3-C4-C5","N1-C5-C4-N3","N1-C2-N3-C4", "H1-C2-N1-C5","H1-C2-N3-C4", "H2-C4-C5-H3", "C4-C5-N1-C7", "C5-C4-N3-C6",
                    "H2-C4-N3-C2", "N1-C5-C4-H2", "N3-C4-C5-H3", "N1-C2-N3-C6", "N3-C2-N1-C7", "H1-C2-N3-C6", "H1-C2-N1-C7",
                    "H2-C4-N3-C6", "H3-C5-N1-C7", "C2-N1-C7-H7","C2-N1-C7-H8","C2-N3-C6-H4","C2-N3-C6-H5","C2-N3-C6-H6",
                    "C4-N3-C6-H4", "C4-N3-C6-H5", "C4-N3-C6-H6", "C5-N1-C7-H8", "C5-N1-C7-H7", "C2-N3-C7-C8", "C2-N1-C7-C8",
                    "C5-N1-C7-C8", "N1-C7-C8-C9", "C8-C9-C10-H13", "C8-C9-C10-H14", "C8-C9-C10-H15", "H11-C9-C10-H13", "H11-C9-C10-H14",
                    "H11-C9-C10-H15", "H12-C9-C10-H13", "H12-C9-C10-H14", "H12-C9-C10-H15", "N1-C7-C8-H9", "N1-C7-C8-H10", "C7-C8-C9-C10",
                    "H7-C7-C8-H9", "H8-C7-C8-H10", "C7-C8-C9-H11", "C7-C8-C9-H12", "H10-C8-C9-C10", "H10-C8-C9-H11", "H10-C8-C9-H12",
                    "H9-C8-C9-C10", "H9-C8-C9-H11", "H9-C8-C9-H12"]
            }
    dihs.update({
            "C2-N3-C4-C5":[1,14.0,2,180,atms["name"].index(dihs['name'][0][:2]),atms["name"].index(dihs['name'][0][3:5]),atms["name"].index(dihs['name'][0][6:8]),atms["name"].index(dihs['name'][0][9:])],
            "N1-C5-C4-N3":[1,14.0,2,180,atms["name"].index(dihs['name'][1][:2]),atms["name"].index(dihs['name'][1][3:5]),atms["name"].index(dihs['name'][1][6:8]),atms["name"].index(dihs['name'][1][9:])],
            "N1-C2-N3-C4":[1,14.0,2,180,atms["name"].index(dihs['name'][2][:2]),atms["name"].index(dihs['name'][2][3:5]),atms["name"].index(dihs['name'][2][6:8]),atms["name"].index(dihs['name'][2][9:])],
            "H1-C2-N1-C5":[2,3.0,2,180,atms["name"].index(dihs['name'][3][:2]),atms["name"].index(dihs['name'][3][3:5]),atms["name"].index(dihs['name'][3][6:8]),atms["name"].index(dihs['name'][3][9:])],
            "H1-C2-N3-C4":[2,3.0,2,180,atms["name"].index(dihs['name'][4][:2]),atms["name"].index(dihs['name'][4][3:5]),atms["name"].index(dihs['name'][4][6:8]),atms["name"].index(dihs['name'][4][9:])],
            "H2-C4-C5-H3":[3,2.0,2,180,atms["name"].index(dihs['name'][5][:2]),atms["name"].index(dihs['name'][5][3:5]),atms["name"].index(dihs['name'][5][6:8]),atms["name"].index(dihs['name'][5][9:])],
            "C4-C5-N1-C7":[4,0.0,1,0,atms["name"].index(dihs['name'][6][:2]),atms["name"].index(dihs['name'][6][3:5]),atms["name"].index(dihs['name'][6][6:8]),atms["name"].index(dihs['name'][6][9:])],
            "C5-C4-N3-C6":[4,0.0,1,0,atms["name"].index(dihs['name'][7][:2]),atms["name"].index(dihs['name'][7][3:5]),atms["name"].index(dihs['name'][7][6:8]),atms["name"].index(dihs['name'][7][9:])],
            "H2-C4-N3-C2":[2,3.0,2,180,atms["name"].index(dihs['name'][8][:2]),atms["name"].index(dihs['name'][8][3:5]),atms["name"].index(dihs['name'][8][6:8]),atms["name"].index(dihs['name'][8][9:])],
            "N1-C5-C4-H2":[2,3.0,2,180,atms["name"].index(dihs['name'][9][:2]),atms["name"].index(dihs['name'][9][3:5]),atms["name"].index(dihs['name'][9][6:8]),atms["name"].index(dihs['name'][9][9:])],
            "N3-C4-C5-H3":[2,3.0,2,180,atms["name"].index(dihs['name'][10][:2]),atms["name"].index(dihs['name'][10][3:5]),atms["name"].index(dihs['name'][10][6:8]),atms["name"].index(dihs['name'][10][9:])],
            "N1-C2-N3-C6":[5,0.0,2,180,atms["name"].index(dihs['name'][11][:2]),atms["name"].index(dihs['name'][11][3:5]),atms["name"].index(dihs['name'][11][6:8]),atms["name"].index(dihs['name'][11][9:])],
            "N3-C2-N1-C7":[5,0.0,2,180,atms["name"].index(dihs['name'][12][:2]),atms["name"].index(dihs['name'][12][3:5]),atms["name"].index(dihs['name'][12][6:8]),atms["name"].index(dihs['name'][12][9:])],
            "H1-C2-N3-C6":[5,0.0,2,180,atms["name"].index(dihs['name'][13][:2]),atms["name"].index(dihs['name'][13][3:5]),atms["name"].index(dihs['name'][13][6:8]),atms["name"].index(dihs['name'][13][9:])],
            "H1-C2-N1-C7":[5,0.0,2,180,atms["name"].index(dihs['name'][14][:2]),atms["name"].index(dihs['name'][14][3:5]),atms["name"].index(dihs['name'][14][6:8]),atms["name"].index(dihs['name'][14][9:])],
            "H2-C4-N3-C6":[5,0.0,2,180,atms["name"].index(dihs['name'][15][:2]),atms["name"].index(dihs['name'][15][3:5]),atms["name"].index(dihs['name'][15][6:8]),atms["name"].index(dihs['name'][15][9:])],
            "H3-C5-N1-C7":[5,0.0,2,180,atms["name"].index(dihs['name'][16][:2]),atms["name"].index(dihs['name'][16][3:5]),atms["name"].index(dihs['name'][16][6:8]),atms["name"].index(dihs['name'][16][9:])],
            "C2-N1-C7-H7":[6,0.195,2,180,atms["name"].index(dihs['name'][17][:2]),atms["name"].index(dihs['name'][17][3:5]),atms["name"].index(dihs['name'][17][6:8]),atms["name"].index(dihs['name'][17][9:])],
            "C2-N1-C7-H8":[6,0.195,2,180,atms["name"].index(dihs['name'][18][:2]),atms["name"].index(dihs['name'][18][3:5]),atms["name"].index(dihs['name'][18][6:8]),atms["name"].index(dihs['name'][18][9:])],
            "C2-N3-C6-H4":[6,0.195,2,180,atms["name"].index(dihs['name'][19][:2]),atms["name"].index(dihs['name'][19][3:5]),atms["name"].index(dihs['name'][19][6:8]),atms["name"].index(dihs['name'][19][9:])],
            "C2-N3-C6-H5":[6,0.195,2,180,atms["name"].index(dihs['name'][20][:2]),atms["name"].index(dihs['name'][20][3:5]),atms["name"].index(dihs['name'][20][6:8]),atms["name"].index(dihs['name'][20][9:])],
            "C2-N3-C6-H6":[6,0.195,2,180,atms["name"].index(dihs['name'][21][:2]),atms["name"].index(dihs['name'][21][3:5]),atms["name"].index(dihs['name'][21][6:8]),atms["name"].index(dihs['name'][21][9:])],
            "C4-N3-C6-H4":[7,0.0,3,0,atms["name"].index(dihs['name'][22][:2]),atms["name"].index(dihs['name'][22][3:5]),atms["name"].index(dihs['name'][22][6:8]),atms["name"].index(dihs['name'][22][9:])],
            "C4-N3-C6-H5":[7,0.0,3,0,atms["name"].index(dihs['name'][23][:2]),atms["name"].index(dihs['name'][23][3:5]),atms["name"].index(dihs['name'][23][6:8]),atms["name"].index(dihs['name'][23][9:])],
            "C4-N3-C6-H6":[7,0.0,3,0,atms["name"].index(dihs['name'][24][:2]),atms["name"].index(dihs['name'][24][3:5]),atms["name"].index(dihs['name'][24][6:8]),atms["name"].index(dihs['name'][24][9:])],
            "C5-N1-C7-H8":[7,0.0,3,0,atms["name"].index(dihs['name'][25][:2]),atms["name"].index(dihs['name'][25][3:5]),atms["name"].index(dihs['name'][25][6:8]),atms["name"].index(dihs['name'][25][9:])],
            "C5-N1-C7-H7":[7,0.0,3,0,atms["name"].index(dihs['name'][26][:2]),atms["name"].index(dihs['name'][26][3:5]),atms["name"].index(dihs['name'][26][6:8]),atms["name"].index(dihs['name'][26][9:])],
            "C2-N3-C7-C8":[8,0.1,3,180,atms["name"].index(dihs['name'][27][:2]),atms["name"].index(dihs['name'][27][3:5]),atms["name"].index(dihs['name'][27][6:8]),atms["name"].index(dihs['name'][27][9:])],
            "C2-N1-C7-C8":[8,0.1,3,180,atms["name"].index(dihs['name'][28][:2]),atms["name"].index(dihs['name'][28][3:5]),atms["name"].index(dihs['name'][28][6:8]),atms["name"].index(dihs['name'][28][9:])],
            "C5-N1-C7-C8":[9,0.2,3,0,atms["name"].index(dihs['name'][29][:2]),atms["name"].index(dihs['name'][29][3:5]),atms["name"].index(dihs['name'][29][6:8]),atms["name"].index(dihs['name'][29][9:])],
            "N1-C7-C8-C9":[7,0.0,3,0,atms["name"].index(dihs['name'][30][:2]),atms["name"].index(dihs['name'][30][3:5]),atms["name"].index(dihs['name'][30][6:8]),atms["name"].index(dihs['name'][30][9:])],
            "C8-C9-C10-H13":[10,0.16,3,0,atms["name"].index(dihs['name'][31][:2]),atms["name"].index(dihs['name'][31][3:5]),atms["name"].index(dihs['name'][31][6:9]),atms["name"].index(dihs['name'][31][10:])],
            "C8-C9-C10-H14":[10,0.16,3,0,atms["name"].index(dihs['name'][32][:2]),atms["name"].index(dihs['name'][32][3:5]),atms["name"].index(dihs['name'][32][6:9]),atms["name"].index(dihs['name'][32][10:])],
            "C8-C9-C10-H15":[10,0.16,3,0,atms["name"].index(dihs['name'][33][:2]),atms["name"].index(dihs['name'][33][3:5]),atms["name"].index(dihs['name'][33][6:9]),atms["name"].index(dihs['name'][33][10:])],
            "H11-C9-C10-H13":[10,0.16,3,0,atms["name"].index(dihs['name'][34][:3]),atms["name"].index(dihs['name'][34][4:6]),atms["name"].index(dihs['name'][34][7:10]),atms["name"].index(dihs['name'][34][11:])],
            "H11-C9-C10-H14":[10,0.16,3,0,atms["name"].index(dihs['name'][35][:3]),atms["name"].index(dihs['name'][35][4:6]),atms["name"].index(dihs['name'][35][7:10]),atms["name"].index(dihs['name'][35][11:])],
            "H11-C9-C10-H15":[10,0.16,3,0,atms["name"].index(dihs['name'][36][:3]),atms["name"].index(dihs['name'][36][4:6]),atms["name"].index(dihs['name'][36][7:10]),atms["name"].index(dihs['name'][36][11:])],
            "H12-C9-C10-H13":[10,0.16,3,0,atms["name"].index(dihs['name'][37][:3]),atms["name"].index(dihs['name'][37][4:6]),atms["name"].index(dihs['name'][37][7:10]),atms["name"].index(dihs['name'][37][11:])],
            "H12-C9-C10-H14":[10,0.16,3,0,atms["name"].index(dihs['name'][38][:3]),atms["name"].index(dihs['name'][38][4:6]),atms["name"].index(dihs['name'][38][7:10]),atms["name"].index(dihs['name'][38][11:])],
            "H12-C9-C10-H15":[10,0.16,3,0,atms["name"].index(dihs['name'][39][:3]),atms["name"].index(dihs['name'][39][4:6]),atms["name"].index(dihs['name'][39][7:10]),atms["name"].index(dihs['name'][39][11:])],
            "N1-C7-C8-H9":[7,0.0,3,0,atms["name"].index(dihs['name'][40][:2]),atms["name"].index(dihs['name'][40][3:5]),atms["name"].index(dihs['name'][40][6:8]),atms["name"].index(dihs['name'][40][9:])],
            "N1-C7-C8-H10":[7,0.0,3,0,atms["name"].index(dihs['name'][41][:2]),atms["name"].index(dihs['name'][41][3:5]),atms["name"].index(dihs['name'][41][6:8]),atms["name"].index(dihs['name'][41][9:])],
            "C7-C8-C9-C10":[11,0.15,1,0,atms["name"].index(dihs['name'][42][:2]),atms["name"].index(dihs['name'][42][3:5]),atms["name"].index(dihs['name'][42][6:8]),atms["name"].index(dihs['name'][42][9:])],
            "H7-C7-C8-H9":[12,0.195,3,0,atms["name"].index(dihs['name'][43][:2]),atms["name"].index(dihs['name'][43][3:5]),atms["name"].index(dihs['name'][43][6:8]),atms["name"].index(dihs['name'][43][9:])],
            "H8-C7-C8-H10":[12,0.195,3,0,atms["name"].index(dihs['name'][44][:2]),atms["name"].index(dihs['name'][44][3:5]),atms["name"].index(dihs['name'][44][6:8]),atms["name"].index(dihs['name'][44][9:])],
            "C7-C8-C9-H11":[12,0.195,3,0,atms["name"].index(dihs['name'][45][:2]),atms["name"].index(dihs['name'][45][3:5]),atms["name"].index(dihs['name'][45][6:8]),atms["name"].index(dihs['name'][45][9:])],
            "C7-C8-C9-H12":[12,0.195,3,0,atms["name"].index(dihs['name'][46][:2]),atms["name"].index(dihs['name'][46][3:5]),atms["name"].index(dihs['name'][46][6:8]),atms["name"].index(dihs['name'][46][9:])],
            "H10-C8-C9-C10":[12,0.195,3,0,atms["name"].index(dihs['name'][47][:3]),atms["name"].index(dihs['name'][47][4:6]),atms["name"].index(dihs['name'][47][7:9]),atms["name"].index(dihs['name'][47][10:])],
            "H10-C8-C9-H11":[12,0.195,3,0,atms["name"].index(dihs['name'][48][:3]),atms["name"].index(dihs['name'][48][4:6]),atms["name"].index(dihs['name'][48][7:9]),atms["name"].index(dihs['name'][48][10:])],
            "H10-C8-C9-H12":[12,0.195,3,0,atms["name"].index(dihs['name'][49][:3]),atms["name"].index(dihs['name'][49][4:6]),atms["name"].index(dihs['name'][49][7:9]),atms["name"].index(dihs['name'][49][10:])],
            "H9-C8-C9-C10":[12,0.195,3,0,atms["name"].index(dihs['name'][50][:2]),atms["name"].index(dihs['name'][50][3:5]),atms["name"].index(dihs['name'][50][6:8]),atms["name"].index(dihs['name'][50][9:])],
            "H9-C8-C9-H11":[12,0.195,3,0,atms["name"].index(dihs['name'][51][:2]),atms["name"].index(dihs['name'][51][3:5]),atms["name"].index(dihs['name'][51][6:8]),atms["name"].index(dihs['name'][51][9:])],
            "H9-C8-C9-H12":[12,0.195,3,0,atms["name"].index(dihs['name'][52][:2]),atms["name"].index(dihs['name'][52][3:5]),atms["name"].index(dihs['name'][52][6:8]),atms["name"].index(dihs['name'][52][9:])]
            })
    imps = {
            'name':["H1-C2-N1-N3","C7-N1-C2-C5","C4-N3-C2-C6","H2-C4-C5-N3","H3-C5-N1-C4"]
            }
    imps.update({
            "H1-C2-N1-N3":[1,0.5,0,atms["name"].index(imps['name'][0][:2]),atms["name"].index(imps['name'][0][3:5]),atms["name"].index(imps['name'][0][6:8]),atms["name"].index(imps['name'][0][9:])],
            "C7-N1-C2-C5":[2,0.6,0,atms["name"].index(imps['name'][1][:2]),atms["name"].index(imps['name'][1][3:5]),atms["name"].index(imps['name'][1][6:8]),atms["name"].index(imps['name'][1][9:])],
            "C4-N3-C2-C6":[2,0.6,0,atms["name"].index(imps['name'][2][:2]),atms["name"].index(imps['name'][2][3:5]),atms["name"].index(imps['name'][2][6:8]),atms["name"].index(imps['name'][2][9:])],
            "H2-C4-C5-N3":[1,0.5,0,atms["name"].index(imps['name'][3][:2]),atms["name"].index(imps['name'][3][3:5]),atms["name"].index(imps['name'][3][6:8]),atms["name"].index(imps['name'][3][9:])],
            "H3-C5-N1-C4":[1,0.5,0,atms["name"].index(imps['name'][4][:2]),atms["name"].index(imps['name'][4][3:5]),atms["name"].index(imps['name'][4][6:8]),atms["name"].index(imps['name'][4][9:])]
            })
    # Write to the packmol file
    pmolfile = 'system.pmol'
    pm = open(pmolfile, 'a')
    pm.write("structure tmp/bmim.xyz\n")
    pm.write("  number %s\n" % num_spec)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
    pm.write("end structure\n")
    pm.write("\n")
    pm.close()
    # Store the number of characteristics
    nchar = [25,25,35,53,5]
    ntyps = [11,14,17,12,2]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
