#!/usr/bin/env python
"""
This is a python program that can be used in tandem with the needblock=.false. 
option in the msd code in order to tabulate the number of blocks.
"""

if __name__ == "__main__":
    import numpy as np
    import sys
    from scipy import stats

    if "-h" in sys.argv:
        print("Usage: python parse_msds.py -f [filestring] -n [nblocks]")
        print("filestring is what precedes the block number")
        print("nblocks is the number of blocks")
        print("Note: you should only use this if you ran the calc_msd code")
        print("      with the needblocks=.false. flag")

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

    t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

    # Get number of columns
    with open(filestring+"_0.dat") as f:
        line = f.readline()
        natms = len(line.split())-1
    time = []
    cm_msd = []
    msd_atm = []
    for i in range(natms):
        msd_atm.append([])
    for i in range(nblocks):
        cmfname = "cm_"+filestring+"_"+str(i)+".dat"
        fname   = filestring+"_"+str(i)+".dat" 
        # read in the cm_msd and time
        time.append(np.genfromtxt(cmfname, usecols = 0, unpack=True))
        cm_msd.append(np.genfromtxt(cmfname, usecols = 1, unpack=True))
        for j in range(natms):
            msd_atm[j].append(np.genfromtxt(fname, usecols = j+1, unpack=True))
    
    cm_msd_av  = np.array(cm_msd).mean(axis=0)
    cm_msd_err = np.array(cm_msd).std(axis=0)*t_val/np.sqrt(nblocks)

    g = open(filestring+"and_err.dat", 'w')
    for i in range(len(time[0])):
        g.write("%.2f " % time[0][i])
        for j in range(natms):
            g.write(" %.6f %.6f" % (np.array(msd_atm[j]).mean(axis=0)[i], np.array(msd_atm[j]).std(axis=0)[i]*t_val/np.sqrt(nblocks)))
        g.write("\n")
    g.close()

    np.savetxt("cm_"+filestring+"and_err.dat", np.c_[time[0], cm_msd_av, cm_msd_err])

    tmpmsd = []
    tmpmsd.append(time[0])
    cmtmpmsd=[]
    cmtmpmsd.append(time[0])
    for j in range(natms):
        for i in range(nblocks):
            tmpmsd.append(msd_atm[j][i])
    for i in range(nblocks):
        cmtmpmsd.append(cm_msd[i])
    np.savetxt(filestring+".dat", np.transpose(tmpmsd ))
    np.savetxt("cm_"+filestring+".dat", np.transpose(cmtmpmsd))


    
    print("All MSDs have been parsed")



       

