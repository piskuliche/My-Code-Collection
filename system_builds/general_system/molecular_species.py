from importlib import import_module as moduleload
def molecule(spectype, num_spec, blength):
    defmolec = moduleload('molecules.'+spectype.lower())
    nchar, ntyps, atms, bnds, angs, dihs, imps = defmolec.define_molec(num_spec, blength)
    # Open Packmol File
    pmolfile = 'system.pmol'
    pm = open(pmolfile, 'a')
    pm.write("structure tmp/"+spectype.lower()+".xyz\n")
    pm.write("  number %s\n" % num_spec)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (blength, blength, blength))
    pm.write("end structure\n")
    pm.write("\n")
    pm.close()
    return nchar, ntyps, atms, bnds, angs, dihs, imps

