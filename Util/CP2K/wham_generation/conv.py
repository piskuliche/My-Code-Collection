#!/usr/bin/env python
"""
This is a python code that is used to go into each directory for which a wham window has been calculated and convert it to units of bohr to match the force constant.
"""
if __name__ == "__main__":
    import numpy as np
    import sys

    if len(sys.argv) == 1:
        print("Usage: python conv.py -basename [name] -outname [outname] -n [nwindows] -skip [nskip]")
        sys.exit("For more information run program with '-h' option")
    if "-h" in sys.argv:
        print("Usage: python conv.py -basename [name] -n [nwidows] -outname [outname] -skip [nskip]")
        print("This code converts the distances in each window file to bohrs from angstroms.")
        print("basename: filename in subdirectory")
        print("outname: filename to output in subdirectory")
        print("nwindows: number of windows (which should go from 0 to nwindows-1")
        print("skip: number of lines to skip at start of file")
        sys.exit("Exiting")
    if "-basename" in sys.argv:
        index = sys.argv.index("-basename")+1
        basename = str(sys.argv[index])
    else:
        sys.exit("basename must be specified")
    if "-n" in sys.argv:
        index = sys.argv.index("-n") + 1
        nwins = int(sys.argv[index])
    else:
        sys.exit("nwindows must be specified")
    if "-outname" in sys.argv:
        index = sys.argv.index("-outname") + 1
        outname = str(sys.argv[index])
    else:
        sys.exit("outname must be specified")
    if "-skip" in sys.argv:
        index = sys.argv.index("-skip") + 1
        skip = int(sys.argv[index])
    else:
        skip = 0
    for i in range(nwins):
        print("Working on bin %d" % i)
        step, col = np.genfromtxt(str(i)+'/'+basename, usecols = (0,1), skip_header=skip, unpack = True)
        col  *= 0.529177 # Converts to bohr from angstrom
        np.savetxt(str(i)+'/'+outname, np.c_[step,col])
    print("There are %d values per window" % len(step))


