#!/usr/bin/env python
#This is a simple code that takes nbins and histograms a common file between all of them.


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    

    if len(sys.argv) == 1:
        print("Usage: python wham_histogram.py -basename [name] -n [nwindows] -unit [unit] -nbin [nbins]")
        sys.exit("For more information run program with '-h' option")
    if "-h" in sys.argv:
        print("Usage: python wham_histogram.py -basename [name] -n [nwidows] -unit [unit] -nbin [nbins]")
        print("This code histograms the bins for the wham code so that you can view overlap")
        print("Outputs a plot of the probability distributions versus the distance in bohr (or specified unit)")
        print("basename: filename in subdirectory")
        print("nwindows: number of windows (which should go from 0 to nwindows-1")
        print("unit: lets you specify the distance unit, should match basename unit")
        print("nbin: number of histogram bins")
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
    if "-unit" in sys.argv:
        index = sys.argv.index("-unit") + 1
        unit = str(sys.argv[index])
    if "-nbin" in sys.argv:
        index = sys.argv.index("-nbin") + 1
        nbins = int(sys.argv[index])
    else:
        nbins = 20

    print("For %s bins, histogramming %s" % (nbins, basename))

    x = []

    for win in range(nwins):
        print("Bin number %s completed" % win)
        x.append([])
        x[win]=np.genfromtxt(str(win)+'/'+basename,unpack=True, usecols=(1),comments='#')
        plt.hist(x[win], nbins, normed=1,histtype='step')

    plt.xlabel('distance (%s)' % unit)
    plt.ylabel('Probability')
    plt.show()
