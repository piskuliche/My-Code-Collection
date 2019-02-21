#!/usr/bin/env python
'''
Python program for reading in a msd and block data and outputting the block undertainty of the trajectory.
The program also fits the mean squared displacements and prints the diffusion coefficient to a file.
'''
import numpy as np
import sys
import math
import scipy.stats as stats
from scipy.optimize import curve_fit

# Read in parameters
if len(sys.argv) != 7:
    print("Usage: fit_msds.py file nblocks molname startskip endskip prepend") 
    exit(0)
filename = str(sys.argv[1])
nblocks  = int(sys.argv[2])
molname  = str(sys.argv[3])
startskip = int(sys.argv[4])
endskip   = int(sys.argv[5])
prepend   = str(sys.argv[6])
t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

# Define the Linear Fitting Function
def linear(x, m, b):
    return m*x + b

# Read in the data
time, msd = np.genfromtxt(filename, usecols = (0,1), unpack=True)
blockmsd = []
for i in range(nblocks):
    j = 2 + i
    blockmsd.append(np.genfromtxt(filename, usecols = j, unpack=True))

# Calculate block uncertainty, output to a file.
bl_std = np.std(blockmsd,axis=0)
bl_err = t_val*bl_std
np.savetxt(molname+'_'+'msd.dat', np.c_[time, msd, bl_err])

# Fit the correlation functions
endskip = len(msd) - endskip
if endskip < 0:
    print("Error: Endskip is too big")
    exit(1)
popt, pcov = curve_fit(linear, time[startskip:endskip], msd[startskip:endskip])
popt_bl = []
pcov_bl = []

fittedout = open('fitmsd_'+molname+'.log','w')
for i in time:
    fittedout.write("%s %s\n" % (time, linear(popt[0], popt[1], i)))
fittedout.close()


for i in range(nblocks):
    tmp_popt, tmp_pcov = curve_fit(linear, time[startskip:endskip], blockmsd[i][startskip:endskip])
    popt_bl.append(tmp_popt[0])
d_tot = popt[0]/6*10**(-4)
d_std = np.std(popt_bl)*t_val/6*10**(-4)

output = open('msd_'+molname+'.log', 'a')
output.write("%s %s %s" % (prepend, d_tot, d_std))
output.close()
