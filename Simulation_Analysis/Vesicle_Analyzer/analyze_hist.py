#!/usr/bin/env python
import numpy as np
import sys

nbins = int(sys.argv[1])

def hist_data(filename,nbins,counts):
    x, hist = np.genfromtxt(filename, usecols=(0,1),dtype=float, unpack=True)
    # Reshape hist by # bin axis
    hist = np.reshape(hist, (-1, nbins))
    singleX = np.reshape(x, (-1, nbins))[0]
    if counts == 0:
        counts = np.sum(hist,axis=1)[0]
    totalhist=np.sum(np.divide(hist.T,counts),axis=1)
    np.savetxt("final_"+filename, np.c_[singleX,totalhist])
    return counts


cnts=hist_data("dr_out_hist_raw.dat", nbins, 0)
hist_data("hdr_out_hist_raw.dat", nbins, cnts)
hist_data("ldr_out_hist_raw.dat", nbins, cnts)

