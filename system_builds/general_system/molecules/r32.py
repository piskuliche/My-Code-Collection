def define_molec(num_spec, blength):
    R32_file = "R32.xyz"
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
        'name':["C1","F1","F2","H1","H2"],
        'atype':[1,2,2,3,3],
        'q':[0.0, 0.0, 0.0, 0.0, 0.0 ],
        'kjeps':[0.0, 0.0, 0.0, 0.0, 0.0],
        'rmin':[0.0, 0.0, 0.0, 0.0, 0.0 ],
        'mass':[12.011, 18.997,18.997,1.008, 1.008],
        }
    bnds = {
            'name':["C1-F1", "C1-F2", "C1-H1","C1-H2"]
            }
    bnds.update({
            "C1-F1":[1, 220.0, 1.400,atms["name"].index(bnds['name'][0][:2]),atms["name"].index(bnds['name'][0][3:])],
            "C1-F2":[1, 220.0, 1.400,atms["name"].index(bnds['name'][1][:2]),atms["name"].index(bnds['name'][1][3:])],
            "C1-H1":[2, 400.0, 1.382,atms["name"].index(bnds['name'][2][:2]),atms["name"].index(bnds['name'][2][3:])], 
            "C1-H2":[2, 400.0, 1.382,atms["name"].index(bnds['name'][3][:2]),atms["name"].index(bnds['name'][3][3:])]
            })
    angs = {
            'name':["H1-C1-H2", "F1-C1-F2", "F1-C1-H1", "F1-C1-H2", "F2-C1-H1", "F2-C1-H2"]
            }
    angs.update({
            "H1-C1-H2":[1,15.0,109.5,atms["name"].index(angs['name'][0][:2]),atms["name"].index(angs['name'][0][3:5]),atms["name"].index(angs['name'][0][6:])],
            "F1-C1-F2":[2,50.0,109.5,atms["name"].index(angs['name'][1][:2]),atms["name"].index(angs['name'][1][3:5]),atms["name"].index(angs['name'][1][6:])],
            "F1-C1-H1":[3,30.0,109.5,atms["name"].index(angs['name'][2][:2]),atms["name"].index(angs['name'][2][3:5]),atms["name"].index(angs['name'][2][6:])],
            "F1-C1-H2":[3,30.0,109.5,atms["name"].index(angs['name'][3][:2]),atms["name"].index(angs['name'][3][3:5]),atms["name"].index(angs['name'][3][6:])],
            "F2-C1-H1":[3,30.0,109.5,atms["name"].index(angs['name'][4][:2]),atms["name"].index(angs['name'][4][3:5]),atms["name"].index(angs['name'][4][6:])],
            "F2-C1-H2":[3,30.0,109.5,atms["name"].index(angs['name'][5][:2]),atms["name"].index(angs['name'][5][3:5]),atms["name"].index(angs['name'][5][6:])]
            })
    dihs = {
            'name':[]
            }
    imps = {
            'name':[]
            }

    # Open Packmol File
    pmolfile = 'IL_system.pmol'
    pm = open(pmolfile, 'a')
    pm.write("structure R32.xyz\n")
    pm.write("  number %s\n" % num_spec)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
    pm.write("end structure\n")
    pm.write("\n")
    pm.close() 

    # Write characteristics to file
    nchar = [5,4,6,0,0]
    ntyps = [3,2,3,0,0]
    return nchar, ntyps, atms, bnds, angs, dihs, imps
