#!/usr/bin/env python
"""
This is a python program aimed at calculating widom insertions by creating molecule files for lammps that can then be inserted
General Format:
    Loop over Configurations:
        read in volume
        Loop over Insertions:
            Attempt Insertion
            Use lammps to calculate insertion energy 
            Calculate exponential factor
        calculate average exponential factor over insertions
        calculate excess chemical potential
        calculate Henry's law Coeff
    print out data
"""

def read_log(logfile,c1,c2):
    e, vol = [], []
    with open(logfile) as log:
        readflag = 0
        for line in log:
            if "Loop" in line:
                readflag = 0
            if readflag == 1:
                e.append(float(line.split()[c1]))
                vol.append(float(line.split()[c2]))
            if "Step" in line:
                readflag = 1
    return e, vol

def read_dump():
    cnf_name = "dump.configs"
    T_CN, X_CN, Y_CN, Z_CN = [], [], [], []
    IX, IY, IZ = [], [], []
    with open(cnf_name, 'r') as cnf:
        readflag = 0
        for line in cnf:
            if "ITEM: ATOMS" in line:
                readflag = 1
                t = []
                x = []
                y = []
                z = []
                ix = []
                iy = []
                iz = []
            if readflag == 1 and "ITEM" not in line:
                t.append(int(line.split()[1]))
                x.append(float(line.split()[2]))
                y.append(float(line.split()[3]))
                z.append(float(line.split()[4]))
                ix.append(float(line.split()[5]))
                iy.append(float(line.split()[6]))
                iz.append(float(line.split()[7]))
            if readflag == 1 and "ITEM" in line and "ATOMS" not in line:
                readflag = 0
                T_CN.append(t)
                X_CN.append(x)
                Y_CN.append(y)
                Z_CN.append(z)
                IX.append(ix)
                IY.append(iy)
                IZ.append(iz)

    return T_CN, X_CN, Y_CN, Z_CN, IX, IY, IZ

def write_dump(step, L, boxshift, t, x, y, z, ix, iy, iz, newr, newt, job):
    natms = len(x)
    na    = len(newt)
    shift = L*boxshift
    xlo, xhi, ylo, yhi, zlo, zhi = 0.0, L, 0.0, L, 0.0, L
    maxt = max(t)
    if step == 0:
        wr = open('dump_traj.widom'+str(job),'w')
        wr.write("ITEM: TIMESTEP\n")
        wr.write("%d\n" % step)
        wr.write("ITEM: NUMBER OF ATOMS\n")
        wr.write("%d\n" % int(natms+na))
        wr.write("ITEM: BOX BOUNDS pp pp pp\n")
        wr.write("%.8f %.8f\n" % (xlo, xhi))
        wr.write("%.8f %.8f\n" % (ylo, yhi))
        wr.write("%.8f %.8f\n" % (zlo, zhi))
        wr.write("ITEM: ATOMS id typ x y z ix iy iz\n")
        for a in range(natms):
            wr.write("%d %d %.6f %.6f %.6f %d %d %d\n" % (a+1, t[a], x[a]+shift, y[a]+shift, z[a]+shift, ix[a], iy[a], iz[a]))
        for a in range(len(newr)):
            wr.write("%d %d %.6f %.6f %.6f %d %d %d\n" % (a+1+natms, newt[a]+maxt, newr[a][0], newr[a][1], newr[a][2], 0, 0, 0))
        wr.close()
    else:
        wr = open('dump_traj.widom'+str(job),'a')
        wr.write("ITEM: TIMESTEP\n")
        wr.write("%d\n" % step)
        wr.write("ITEM: NUMBER OF ATOMS\n")
        wr.write("%d\n" % int(natms+na))
        wr.write("ITEM: BOX BOUNDS pp pp pp\n")
        wr.write("%.8f %.8f\n" % (xlo, xhi))
        wr.write("%.8f %.8f\n" % (ylo, yhi))
        wr.write("%.8f %.8f\n" % (zlo, zhi))
        wr.write("ITEM: ATOMS id typ x y z ix iy iz\n")
        for a in range(natms):
            wr.write("%d %d %.6f %.6f %.6f %d %d %d\n" % (a+1, t[a], x[a]+shift, y[a]+shift, z[a]+shift,ix[a], iy[a], iz[a]))
        for a in range(len(newr)):
            wr.write("%d %d %.6f %.6f %.6f %d %d %d\n" % (a+1+natms, newt[a]+maxt, newr[a][0], newr[a][1], newr[a][2], 0, 0, 0))
        wr.close()
        

if __name__ == "__main__":
    import numpy as np
    import gen_ins as gen
    import conv_connectivity as cc
    import sys, json
    # Read the NML file

    if len(sys.argv) > 1 and sys.argv[1] == "-h":
        print("Need to read an input file as python widom_backbone.py < inputfile")
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
        njob     = int(sys.argv[index+1])
    else:
        print("No farm info provided")
        print("If not farming and only have single job run with -farm 1 1")
        sys.exit()
        
    try:
        with open(inpfile) as f:
            nml = json.load(f)
    except json.JSONDecodeError:
        print('Exiting on Invalid JSON format')
        print('run program with "-h" flag to see input file options')
        sys.exit()

    defaults = {"nfile":"traj.xyz", "nconfigs":10, "ninsert":100, "bshift":0.0, "log":"log.production", "pecolumn":3, "volcolumn":5, "conn_file":"test.connect","sconfig":0, "econfig":10,"temp":298.15, "logins":"log.rerun", "num_blocks":5,"num_solvent":343}

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


    # Read the energy of each configuration
    e_configs, v_configs = read_log(logfile, ecol, volcol)
    assert len(e_configs) == num_configs, "There are more log steps (%d) than configs defined (%d)" % (len(e_configs), num_configs)
    # Read in the x y z coordinates of each configuration
    TYP, X,Y,Z, IX, IY, IZ = read_dump()
    print("There are %d Configurations" % len(X))
    # Read the connectivity file
    na, M, ntyp, header, coords, footer = cc.read_connectivity(connectfile)
    print("There are %d atoms in new molecule" % na)
    tmpx, tmpy, tmpz = cc.read_xyz(connectfile)
    tmpr = [np.transpose([tmpx,tmpy,tmpz])]
    # Calculate COM and subtract
    ctmpr = gen.calculate_com(tmpr, M)
    for i in range(len(tmpr)):
        for j in range(3):
            tmpr[i][j] -= ctmpr[j]
    # initialize random generator
    gen.init_random()
    # Loop to make new trajectory file
    count = 0
    if farmjobs > 1:
        npfarm = (end_config - startconfig)/farmjobs
        offset = start_config
        start_config = njob*npfarm + offset
        end_config = (njob+1)*npfarm + offset
    for c in range(start_config, end_config):
        # Loop over insertions
        L = v_configs[c]**(1/3.)
        print("Loop step %d of %d" % (c+1, npfarm+star_config))
        for i in range(num_insert):
            if na > 1:
                r = gen.rand_rot(tmpr)
            else:
                r = tmpr
            r_translated, com_r = gen.rand_trans(r, L)
            write_dump(count, L,boxshift, TYP[c], X[c], Y[c], Z[c],IX[c], IY[c], IZ[c], r_translated, ntyp, njob)
            count += 1
        



