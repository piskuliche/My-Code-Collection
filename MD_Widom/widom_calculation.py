#!/usr/bin/env python
"""This is a code that takes the output of the lammps rerun logfile and calculates the excess chemical potential and the Henry's law coefficient"""

def calc_exp_factor(beta, e_ins):
    """
    Calculates averaged insertion energy boltzmann factor.
    Input:  beta=1/(kb*T)
            e_ins - array with insertion energies for particular config
    Output: efact - average value of the boltzmann factor of insertion energies  
    """
    efact = 0
    for nu in e_ins:
        efact += np.exp(-beta * nu)
    efact /= float(len(e_ins))
    return efact

def calc_muex_H(beta, numerator, denominator, nblocks, nsolv):
    """
    Calculates the excess chemical potential and H-Law Constant
    Input: beta = 1/(kb*T)
           numerator - array that holds all the boltz factors weighted by volume
           denominator - average volume of range
    Output: muex, muex_error, H, H_error
    """
    # Uncertainty
    t_val = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)
    # Calc ratio of averages
    ratio = np.average(numerator)/np.average(denominator)
    rho   = nsolv/np.average(denominator)
    # Calc excess chem pot, H
    muex  = -1/beta*np.log(ratio)
    H     = rho/beta*ratio
    print muex, H

    # Calculate the number of configs per block
    nvals = len(numerator)
    valspblock = nvals/float(nblocks)
    assert nvals % valspblock == 0, "Number of configs (%d) not evenly divisible by number of blocks(%d)!" % (nvals, nblocks)
    
    # Block Averaging
    bl_muex, bl_H = [],[]
    for block in range(nblocks):
        s = int(block*valspblock)
        e = int((block+1)*valspblock)
        ratio_bl = np.average(numerator[s:e])/np.average(denominator[s:e])
        rho = nsolv/np.average(denominator[s:e])
        bl_muex.append(-1/beta*np.log(ratio_bl))
        bl_H.append(rho/beta*ratio_bl)
    print(bl_H) 
    muex_error, H_error = np.std(bl_muex)*t_val/np.sqrt(nblocks-1), np.std(bl_H)*t_val/np.sqrt(nblocks-1)
    return muex, muex_error, H, H_error

 

if __name__ == "__main__":
    import numpy as np
    import gen_ins as gen
    import conv_connectivity as cc
    import widom_backbone as wb
    import sys, json
    from scipy import stats
    # Read the NML file

    if len(sys.argv) > 1 and sys.argv[1] == "-h":
        print("Need to read an input file as python widom_calculate.py -in inputfile -farm njobs")
        print("Possible input options are:")
        print("     nfile - [string]  file that holds the lammps style dumps (id type x y z)")
        print("     nconfigs   - [integer] number of configs in said file")
        print("     ninsert    - [integer] number of attempted insertions per config")
        print("     bshift     - [float]   shifts box by a multiplicative factor of the box length (i.e. 0.5 adds 0.5 L to each coord)")
        print("     logrun     - [string]  log file by lammps - must only include dumps from the configs")
        print("     pecolumn   - [integer] column value of the potential energy in log file")
        print("     volcolumn  - [integer] column value of the volume in the log file")
        print("     conn_file  - [string]  lammps connectivity file for inserted molecule")
        print("     sconfig    - [integer]  starting configuration number")
        print("     econfig    - [integer]  ending configuration number")
        print("     temp       - [float]    simulation temperature")
        print("     logins     - [string]   log file from lammps, most only include dumps from the insertions")
        print("     num_blocks    - [integer]  number of blocks for block averaging")
        print("     num_solvent   - [integer]  number of solvent molecules")
        sys.exit()
    if "-in" in sys.argv:
        index = sys.argv.index("-in")+1
        inpfile = str(sys.argv[index])
    else:
        print("No input file provided")
        sys.exit()
    if "-farm" in sys.argv:
        index    = sys.argv.index("-farm")+1
        farmjobs = int(sys.argv[index])
    else:
        print("No farm info provided")
        print("If not farming and only have single job run with -farm 1 ")
        sys.exit()

    try:
        with open(inpfile) as f:
            nml = json.load(f)
    except json.JSONDecodeError:
        print('Exiting on Invalid JSON format')
        print('run program with "-h" flag to see input file options')
        sys.exit()

    defaults = {"nfile":"traj.xyz", "nconfigs":10, "ninsert":100, "bshift":0.0, "log":"log.production", "pecolumn":3, "volcolumn":5, "conn_file":"test.connect","sconfig":0, "econfig":10, "temp":298.15, "logins":"log.rerun", "num_blocks":5, "num_solvent":343}

    configfile      = nml["nfile"]              if "nfile"              in nml else defaults["nfile"]
    num_configs     = nml["nconfigs"]           if "nconfigs"           in nml else defaults["nconfigs"]
    num_insert      = nml["ninsert"]            if "ninsert"            in nml else defaults["ninsert"]
    boxshift        = nml["bshift"]             if "bshift"             in nml else defaults["bshift"]
    logfile         = nml["logrun"]             if "logrun"             in nml else defaults["logrun"]
    ecol            = nml["pecolumn"]           if "pecolumn"           in nml else defaults["pecolumn"]
    volcol          = nml["volcolumn"]          if "volcolumn"          in nml else defaults["volcolumn"]
    connectfile     = nml["conn_file"]          if "conn_file"          in nml else defaults["conn_file"]
    start_config    = nml["sconfig"]            if "sconfig"            in nml else defaults["sconfig"]
    end_config      = nml["econfig"]            if "econfig"            in nml else defaults["econfig"]
    temperature     = nml["temp"]               if "temp"               in nml else defaults["temp"]
    ins_logfile     = nml["logins"]             if "logins"             in nml else defaults["logins"]
    nblocks         = nml["num_blocks"]         if "num_blocks"         in nml else defaults["num_blocks"]
    nsolvent        = nml["num_solvent"]        if "num_sovlent"        in nml else defaults["num_solvent"]
    
    # Define constants
    kb = 0.0019872041 # kcal/(mol K)
    beta = 1.0/(kb*temperature) # mol/kcal
    # Read in the insertion energies
    if farmjobs == 1:
        pe, vol = wb.read_log(ins_logfile+str(1), ecol, volcol)
    else:
        pe, vol = [], []
        for i in range(farmjobs):
            tmpp, tmpv = wb.read_log(ins_logfile+str(i), ecol, volcol)
            for val in tmpp:
                pe.append(val)
            for val in tmpv:
                vol.append(val)
    # Read in the config energies
    pe_conf, vol = wb.read_log(logfile, ecol, volcol)
    

    boltz_fact = []
    num, denom = [], []
    for c in range(start_config, end_config):
        """This calculates the ratio <V <e^-(beta*nu)>>/<V>"""
        s = c*num_insert
        e = (c+1)*num_insert
        print s, e
        # Calculates <e^-(beta*nu)> over insertions in config
        boltz_fact.append(calc_exp_factor(beta, np.array(pe[s:e])-pe_conf[c]))
        num.append(vol[c]*boltz_fact[c])
        denom.append(vol[c])
    ex_chempot, err_ex_chempot, henry, err_henry = calc_muex_H(beta, num, denom, nblocks, nsolvent)
    print(ex_chempot, err_ex_chempot, henry, err_henry)




