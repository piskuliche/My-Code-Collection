#!/usr/bin/env python
"""
This is a code that is meant to plot, and manipulate the gofr's output by the fortran code.
"""
def read_gofr(t_val,selec1, selec2, nblocks):
    # Read in Data
    fstring="pairdist_"+selec1+"_"+selec2+".dat"
    dist, gofr = np.genfromtxt(fstring, usecols=(0,1), unpack=True)

    # Calculates the RDF
    bl_gofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_gofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))

    err_gofr = np.std(bl_gofr,axis=0)*t_val
    return dist, gofr, bl_gofr, err_gofr

def read_egofr(t_val,selec1, selec2, nblocks):
    # Read in Data
    fstring="epairdist_"+selec1+"_"+selec2+".dat"
    dist, egofr = np.genfromtxt(fstring, usecols=(0,1), unpack=True)

    # Calculates the RDF
    """
    bl_egofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_egofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))

    err_egofr = np.std(bl_egofr,axis=0)*t_val
    """
    #return egofr, bl_egofr, err_egofr
    return egofr

def calc_nofr(gofr, dist, bl_gofr, L, N):
    nofr= integrate.cumtrapz(gofr*dist**2.*N/L**3.*4*np.pi, dist, initial=0)
    bl_nofr=[] 
    for i in range(1,nblocks+1):
        bl_nofr.append(integrate.cumtrapz(bl_gofr[i-1]*dist**2.*N/L**3.*4*np.pi, dist, initial=0))
    err_nofr = np.std(bl_nofr,axis=0)*t_val 
    return nofr, err_nofr

def peak_find(dist, gofr, nofr, thresh):
    # Maxima
    mlocs, dict =  find_peaks(gofr, threshold=thresh, height=1)
    maxvals = []
    dmax = []
    for i in mlocs:
        print('Maximum value at %4.5f %4.5f with coordination number %4.5f' %(dist[i],gofr[i],nofr[i]))
        maxvals.append(gofr[i])
        dmax.append(dist[i])
    minvals = []
    dmin = []
    for i in range(len(mlocs)-1):
        loc = np.where(gofr == min(gofr[mlocs[i]:mlocs[i+1]]))[0]
        minvals.append(gofr[loc])
        dmin.append(dist[loc])
        print('Minimum value at %4.5f %4.5f with coordination number %4.5f' %(dist[loc],gofr[loc],nofr[loc]))
    return maxvals, dmax, minvals, dmin

def calc_PMF(gofr, T):
    PMF = []
    pdist=[]
    c=0
    for i in range(len(gofr)):
        if gofr[i] > 0:
            PMF.append(-kb*T*np.log(gofr[i]))
        else:
            c += 1
    return PMF, c

def calc_PMF_deriv(PMF, gofr, egofr, T):
    dPMF = []
    for i in range(len(PMF)):
        dPMF.append(kb*T*(egofr[i]/gofr[i]-PMF[i]))
    return dPMF

def pred_PMF(dist, PMF, dPMF, T, Tp):
    pred = []
    for i in range(len(PMF)):
        pred.append(PMF[i]+dPMF[i]*(1/(kb*Tp)-1/(kb*T)))
    np.savetxt('pmf_pred.dat', np.c_[dist,pred])
    return pred


    


import numpy as np
import argparse, math
import matplotlib.pyplot as plt
from scipy import stats,integrate
from scipy.signal import find_peaks

kb=0.0019872041 #kcal/mol


# Command Line Input
parser = argparse.ArgumentParser()
parser.add_argument('-p', default=0, help='0 Plots in matplotlib')
parser.add_argument('-selec1',default=1, help='selection1')
parser.add_argument('-selec2', default=1, help='selection2')
parser.add_argument('-nblocks',default=5, help='nblocks')
parser.add_argument('-L', default=21.752, help='box side length')
parser.add_argument('-N', default=343, help='Number of atoms')
parser.add_argument('-T', default=298.15, help='Temperature')
parser.add_argument('-Tpred', default=240, help='Temperature to predict')
parser.add_argument('-thresh', help='threshold for peak finding')
parser.add_argument('-ew', default=0, help='weighted e off [0] and on [1]')
args = parser.parse_args()
selec1  = str(args.selec1)
selec2  = str(args.selec2)
nblocks = int(args.nblocks)
L       = float(args.L)
N       = int(args.N)
T       = float(args.T)
Tp      = float(args.Tpred)
ew      = int(args.ew)
pltflag = int(args.p)
if args.thresh is not None:
    thresh=float(args.thresh)
else:
    thresh=0.001

t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)


# Read GofR
dist, gofr, bl_gofr, err_gofr = read_gofr(t_val,selec1, selec2, nblocks)

# Read dHGofR
if ew == 1:
    egofr = read_egofr(t_val, selec1, selec2, nblocks)

# Calculates the Coordination Number
nofr, bl_nofr = calc_nofr(gofr, dist, bl_gofr, L, N)

# Finds the local minima and maxima
maxvals, dmax, minvals, dmin = peak_find(dist, gofr, nofr,thresh)

# Calculates the PMF

PMF, c = calc_PMF(gofr, T)
np.savetxt("pmf_"+str(selec1)+"_"+str(selec2)+".dat", np.c_[dist[c:], PMF])

# Calculates the Derivative of the PMF
if ew == 1:
    dPMF=calc_PMF_deriv(PMF, gofr[c:], egofr[c:], T)
    np.savetxt("dpmf_"+str(selec1)+"_"+str(selec2)+".dat", np.c_[dist[c:], dPMF])

# Makes PMF Prediction
if ew == 1:
    pPMF = pred_PMF(dist[c:], PMF, dPMF, T, Tp)

if pltflag == 1:
    plt.plot(dmax,maxvals,'bs')
    plt.plot(dmin,minvals,'rs')
    plt.plot(dist,gofr,'-b')
    plt.show()
    if ew == 1:
        plt.plot(dist,egofr)
        plt.show()
    plt.plot(dist,nofr)
    plt.show()
    plt.plot(dist[c:],PMF)
    if ew == 1:
        plt.plot(dist[c:],dPMF)
    plt.show()
    



