from importlib import import_module as moduleload
def molecule(spectype, num_spec, blength):
    defmolec = moduleload('molecules.'+spectype.lower())
    nchar, ntyps, atms, bnds, angs, dihs, imps = defmolec.define_molec(num_spec, blength)
    return nchar, ntyps, atms, bnds, angs, dihs, imps

