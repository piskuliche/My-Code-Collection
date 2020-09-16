from importlib import import_module as moduleload
def molecule(spectype, num_spec, blength):
    defmolec = moduleload('molecules.'+spectype.lower())
    nchar, ntyps, atms, bnds, angs, dihs, imps = defmolec.define_molec(num_spec, blength)
    # Open Packmol File
    pmolfile = 'system.pmol'
    pm = open(pmolfile, 'a')
    pm.write("structure tmp/"+spectype.lower()+".xyz\n")
    pm.write("  number %s\n" % num_spec)
    #pm.write("  constrain_rotation x 0. 0.\n")
    #pm.write("  constrain_rotation y 0. 0.\n")
    #pm.write("  constrain_rotation z 0. 0.\n")
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength-2., blength-2., blength-2.))
    pm.write("end structure\n")
    pm.write("\n")
    pm.close()
    return nchar, ntyps, atms, bnds, angs, dihs, imps

