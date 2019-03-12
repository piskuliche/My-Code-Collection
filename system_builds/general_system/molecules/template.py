def define_molec(num_spec, blength):
    mol_file = ".xyz"
        mol = open(mol_file, 'w')
        # XYZ Header
        mol.write("\n")
        mol.write("\n")
        #  LABEL XYZ 
        mol.close()
        atms = {
            'name':[],
            'atype':[],
            'q':[],
            'kjeps':[],
            'rmin':[],
            'mass':[],
            }
        bnds = {
                'name':[]
                }
        bnds.update({
                })
        angs = {
                'name':[]
                }
        angs.update({
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
        pm.write("structure .xyz\n")
        pm.write("  number %s\n" % num_spec)
        pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
        pm.write("end structure\n")
        pm.write("\n")
        pm.close()

        # Write characteristics to file
        nchar = [0,0,0,0,0]
        ntyps = [0,0,0,0,0]

    return nchar, ntyps, atms, bnds, angs, dihs, imps
