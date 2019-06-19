#!/usr/bin/env python
"""
This is a code that is meant to plot, and manipulate the gofr's output by the fortran code.
"""
import numpy as np
import argparse, math
import matplotlib.pyplot as plt
from scipy import stats,integrate



# Command Line Input
parser = argparse.ArgumentParser()
parser.add_argument('-p', default=0, help='0 Plots in matplotlib')
parser.add_argument('-selec1', help='selection1')
parser.add_argument('-selec2', help='selection2')
parser.add_argument('-nblocks', help='nblocks')
args = parser.parse_args()
selec1  = str(args.selec1)
selec2  = str(args.selec2)
nblocks = int(args.nblocks)
pltflag = int(args.p)

t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)


# Read in Data
fstring="pairdist_"+selec1+"_"+selec2+".dat"

dist, gofr = np.genfromtxt(fstring, usecols=(0,1), unpack=True)

# Calculates the RDF
bl_gofr = []
for i in range(1,nblocks+1):
    bstring=("bl_%d_"%i)
    bl_gofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))

err_gofr = np.std(bl_gofr,axis=0)*t_val

# Calculates the Coordination Number
nofr= integrate.cumtrapz(gofr*dist**2.*343/21.752**3.*4*np.pi, dist, initial=0)
bl_nofr=[]
for i in range(1,nblocks+1):
    bl_nofr.append(integrate.cumtrapz(bl_gofr[i-1]*dist**2.*343/21.752**3.*4*np.pi, dist, initial=0))
err_nofr = np.std(bl_nofr,axis=0)*t_val

# Calculates the PMF

if pltflag == 1:
    plt.errorbar(dist,gofr,yerr=err_gofr)
    plt.show()
    plt.errorbar(dist,nofr, yerr=err_nofr)
    plt.show()
    



