def define_molec(num_spec, blength):
    mol_file = "tmp/.xyz"
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

        # Write characteristics to file
        nchar = [0,0,0,0,0]
        ntyps = [0,0,0,0,0]

    return nchar, ntyps, atms, bnds, angs, dihs, imps
