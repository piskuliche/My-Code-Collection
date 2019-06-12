#!/usr/bin/env python
"""
This is a python program that can be used in tandem with the needblock=.false. 
option in the msd code in order to tabulate the number of blocks.
"""

if __name__ == "__main__":
    import numpy as np
    import sys
    from scipy import stats
    from fit_msds import D_from_MSD
    

    usage="Usage: python parse_msds.py -f [filestring] -n [nblocks] -startskip [startskip] -endskip [endskip] -molname [molname] -label [label]"
    if len(sys.argv) == 1:
        print("%s" % usage)
    if "-h" in sys.argv:
        print("%s" % usage)
        print("filestring is what precedes the block number")
        print("nblocks is the number of blocks")
        print("startskip is the number of data points to skip at start in fitting")
        print("endskip is the number of data points to skip at end in fitting")
        print("molname is the molecule name")
        print("label is a label used in the diffusion output file")
        print("Note: you should only use this if you ran the calc_msd code")
        print("      with the needblocks=.false. flag")
        sys.exit()

    if "-f" in sys.argv:
        index = sys.argv.index("-f")+1
        filestring = str(sys.argv[index])
    else:
        print("-f option not included")
        print("run with -h flag for options")
        sys.exit()
    if "-n" in sys.argv:
        index = sys.argv.index("-n")+1
        nblocks = int(sys.argv[index])
    else:
        print("-n option not included")
        print("run with -h flag for options")
        sys.exit()

    if "-startskip" in sys.argv:
        index = sys.argv.index("-startskip")+1
        startskip = int(sys.argv[index])
    else:
        startskip = 0
        print("Defaulting to startskip of 0")

    if "-endskip" in sys.argv:
        index = sys.argv.index("-endskip")+1
        endskip = int(sys.argv[index])
    else:
        endskip = 0
        print("Defaulting to endskip of 0")
    
    if "-molname" in sys.argv:
        index = sys.argv.index("-molname")+1
        molname = str(sys.argv[index])
    else:
        molname = "molecule"
        print("Defaulting to molname of 'molecule'")
    if "-label" in sys.argv:
        index = sys.argv.index("-label")+1
        label = str(sys.argv[index])
    else:
        label = molname
        print("Defaulting to label of %s" % molname)

    t_val=stats.t.ppf(0.975,nblocks-1)

    # Get number of columns
    with open(filestring+"_0.dat") as f:
        line = f.readline()
        natms = len(line.split())-1
    time = []
    cm_msd = []
    msd_atm = []

    # Read in the msds for each atom and the center of mass
    for i in range(natms):
        msd_atm.append([])
    for i in range(nblocks):
        cmfname = "cm_"+filestring+"_"+str(i)+".dat"
        fname   = filestring+"_"+str(i)+".dat" 
        # read in the cm_msd and time
        time.append(np.genfromtxt(cmfname, usecols = 0, unpack=True))
        cm_msd.append(np.genfromtxt(cmfname, usecols = 1, unpack=True))
        # read in the atom msds
        for j in range(natms):
            msd_atm[j].append(np.genfromtxt(fname, usecols = j+1, unpack=True))
    # Calculate the average cm msd
    cm_msd_av  = np.array(cm_msd).mean(axis=0)
    cm_msd_err = np.array(cm_msd).std(axis=0)*t_val/np.sqrt(nblocks)
    
    # Set endskip based on length of msd
    endskip = len(cm_msd_av) - endskip

    # Print the per atom msds
    g = open(filestring+"and_err.dat", 'w')
    for i in range(len(time[0])):
        g.write("%.2f " % time[0][i])
        for j in range(natms):
            g.write(" %.6f %.6f" % (np.array(msd_atm[j]).mean(axis=0)[i], np.array(msd_atm[j]).std(axis=0)[i]*t_val/np.sqrt(nblocks)))
        g.write("\n")
    g.close()

    # Print the CM MSD with error
    np.savetxt("cm_"+filestring+"and_err.dat", np.c_[time[0], cm_msd_av, cm_msd_err])

    """ This section used to print out a huge file - kinda worthless"""
    #tmpmsd = []
    #tmpmsd.append(time[0])
    #cmtmpmsd=[]
    #cmtmpmsd.append(time[0])
    #for j in range(natms):
    #    for i in range(nblocks):
    #        tmpmsd.append(msd_atm[j][i])
    #for i in range(nblocks):
    #    cmtmpmsd.append(cm_msd[i])
    #np.savetxt(filestring+".dat", np.transpose(tmpmsd ))
    #np.savetxt("cm_"+filestring+".dat", np.transpose(cmtmpmsd))

    # Calculates the Diffusion coefficients for the COM
    d_cm, d_std = D_from_MSD(time[0], cm_msd_av, cm_msd, startskip, endskip, nblocks, molname, label, t_val)
    print("Diffusion %s: %2.4e +/- %2.4e" % (label,d_cm, d_std))

    
    print("All MSDs have been parsed")



       

