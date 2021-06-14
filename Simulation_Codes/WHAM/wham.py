#!/usr/bin/env python
import numpy as np
from numpy import inf,nan
import argparse
import matplotlib.pyplot as plt

"""
This is a simple Weighted Histogram Analysis Code that can be used to calculate the Potential of Mean Force.
Please read the readme file for more specific notes about its operation.

Copyright 2021, Zeke Piskulich, University of Kansas
For questions, concerns, email piskulichz@gmail.com
"""

# Constants
kb=0.0019872041
T=298.15
kbT=kb*T

def calc_bias(x,k):
    # This function calculates the bias potential value 
    return 0.5*k*x**2.

def calc_p(F):
    # This function calculates the unweighted probability distribution
    numerator=np.sum(ni,axis=0)
    denominator=np.sum(np.multiply(cnt,np.exp(np.divide(-np.subtract(U.T,F-np.min(F)),kbT))),axis=1)
    return numerator/denominator

def calc_f(P):
    # This calculates the weight of each window (and then takes the log to get the free energy)
    F=-kbT*np.log(np.sum(np.multiply(P,np.exp(-U/kbT)),axis=1))
    F[F==inf] = 0.0
    return F




def wham_iteration(F):
    P=calc_p(F)
    F=calc_f(P)
    # Calculates and writes out the error
    dF = np.sum(np.abs(np.subtract(F,F_old)))
    print("Error: %s " % dF)
    # Checks if convered, if it is, writes files and exits calculation.
    isconverged=False
    if (dF < tolerance):
        isconverged=True
        np.savetxt("wham_pmf.out", np.c_[xvals,-kbT*np.log(P)-np.min(-kbT*np.log(P)),P])
        np.savetxt("wham_F.out", np.c_[F-np.min(F)])
    return F, isconverged

# START

#Command Line Arguments
parser=argparse.ArgumentParser()
parser.add_argument('-Nw', default=30,type=int, help='Number of Windows')
parser.add_argument('-rlow', default=1.5,type=float, help='Low cutoff')
parser.add_argument('-rhigh', default=8,type=float, help='High cutoff')
parser.add_argument('-nbin', default=100, type=int, help='Histogram bins')
parser.add_argument('-k', default=11.0, type=float, help='Force constant')
parser.add_argument('-plot', default=0, type=int, help='If 1, plots histograms in matplotlib')
parser.add_argument('-subfile', default="lif.distance", type=str, help="File name for colvar in subdirectory")
parser.add_argument('-unit', default=1, type=int, help="[0] Angstrom [1] bohr (output always angstrom, kcal/mol)")
args = parser.parse_args()

Nwindows=args.Nw
rlow=args.rlow
rhi=args.rhigh
nbins=args.nbin
k=args.k
unit=args.unit
subfile=args.subfile
shouldplot=args.plot

# Set Units
BohrToAng=0.529177
convdist=1.0
convk=1.0
if unit==1:
    print("Converting r,k to units of angstroms, kcal/mol/ang^2 from bohr, hartree/bohr^2")
    convdist=BohrToAng
    convk=627.509/(BohrToAng**2.)

# Does unit conversions
k=k*convk
rlow=rlow*convdist
rhi=rhi*convdist
print("New low: %s New hi: %s" % (rlow, rhi))
print("new k: %s" % k)

# Sets up Arrays
U=[]
ni=[]
cnt=[]
center=[]
F_old=np.zeros(Nwindows)
xvals=np.linspace(rlow,rhi,num=nbins)
# Reads in bias potential locations
xc=np.genfromtxt("wham_metadata.info",usecols=1,unpack=True)
# Sets up the Histograms
for window in range(Nwindows):
    # Reads Colvar Data
    data=np.genfromtxt(str(window)+"/"+subfile, usecols=1,unpack=True)
    data=data*convdist # Converts data
    # Histograms data
    hist,bins=np.histogram(data,bins=nbins,range=(rlow,rhi),density=False)
    # Plotting Stuff
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    ni.append(hist) # Array that stores counts for each x, window
    cnt.append(np.sum(hist)) # Total number of counts for each window
    dx=np.subtract(xc[window],xvals) # difference between the bias location and each x
    U.append(calc_bias(dx,k)) # Array of the bias potential (2d array of shape dx, nwindow)

U=np.array(U)

# Plots the histograms using matplotlib
if shouldplot==1: plt.show()

# Basic Settings
isconverged = False
tolerance = 1e-8
maxiter=10000
iteration=0
#does the wham iterations
while ( isconverged == False):
    F_old,isconverged=wham_iteration(F_old)
    iteration+=1
    if iteration > maxiter:
        exit("Too many iterations")
    
    

     



