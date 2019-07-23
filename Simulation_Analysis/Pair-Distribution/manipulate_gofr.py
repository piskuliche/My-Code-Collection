#!/usr/bin/env python
"""
This is a code that is meant to plot, and manipulate the gofr's output by the fortran code.
"""

def read_gofr(t_val,selec1, selec2, nblocks):
    """
    This function reads in the gofr data (two column distance gofr formatted)
    from a file (i.e. pairdist_1_1.dat)
    and then the block gofrs from the block flies (two column dist gofr formt)
    from the files (i.e. bl_1_pairdist_1_1.dat)
    It then calculates the error in the RDF from the blocks, according to a 95% confidence
    interval associated with Student's t-distribution
    Writes to a file (i.e RDF_1_1_298.dat)
    """
    # Reads in RDF
    fstring="pairdist_"+selec1+"_"+selec2+".dat"
    dist, gofr = np.genfromtxt(fstring, usecols=(0,1), unpack=True)

    # Reads in Blcok RDF
    bl_gofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_gofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))
    # Calculates the confidence interval (note t_val includes a factor of 1/sqrt(nblocks) implicitly)
    err_gofr = np.std(bl_gofr,axis=0)*t_val
    np.savetxt('RDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,gofr,err_gofr])
    return dist, gofr, bl_gofr, err_gofr

def read_egofr(t_val,selec1, selec2, nblocks):
    """
    This function reads in the gofr data (two column distance gofr formatted)
    from a file (i.e. epairdist_1_1.dat)
    and then the block gofrs from the block flies (two column dist gofr formt)
    from the files (i.e. bl_1_epairdist_1_1.dat)
    It then calculates the error in the RDF from the blocks, according to a 95% confidence
    interval associated with Student's t-distribution
    Writes to a file (i.e eRDF_1_1_298.dat)
    """
    # Read in eRDF
    fstring="epairdist_"+selec1+"_"+selec2+".dat"
    dist, egofr = np.genfromtxt(fstring, usecols=(0,1), unpack=True)

    # Reads in block eRDF
    bl_egofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_egofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))

    err_egofr = np.std(bl_egofr,axis=0)*t_val
    np.savetxt('eRDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,egofr,err_egofr])
    return egofr, bl_egofr, err_egofr

def pred_gofr(gofr, egofr, bl_gofr, bl_egofr, T, Tp):
    """
    This function calculates the first order taylor series prediction based on the gofr and the derivative
    Note: egofr = <dH delta(r-r')> = - d(RDF)/d(beta)
    Thus the TS expansion is:
    RDF(beta) = RDF(beta_o) + d(RDF)/d(beta)|beta_o (beta-beta_o)
              = RDF(beta_o) - <dH delta(r-r')>|beta_o (beta-beta_o)
    This prediction is saved to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    """
    pred = []
    taylorfactor=1/(kb*Tp) - 1/(kb*T) # (beta-beta_o)
    # Makes the prediction for each point in the gofr
    for i in range(len(gofr)):
        pred.append(gofr[i]-egofr[i]*taylorfactor)
    # Block averaging
    bl_pred = []
    for b in range(nblocks):
        tmp_pred = []
        for i in range(len(bl_gofr[b])):
            tmp_pred.append(bl_gofr[b][i]-bl_egofr[b][i]*taylorfactor)
        bl_pred.append(tmp_pred)
    err_pred = np.std(bl_pred, axis=0)*t_val
    # Saves the result to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    if graphmovie == 0:
        np.savetxt('gofr_pred_'+str(selec1)+"_"+str(selec2)+"_"+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_pred])
    else:
        np.savetxt('gofr_frame_'+str(int(Tp))+'.dat', np.c_[dist,pred,err_pred])
    return pred


def calc_nofr(gofr, dist, bl_gofr, L, N):
    """
    This piece of code calculates the coordination number.
    This is achieved via a trapezoid rule integration
    n = (4*pi*N/V)∫<delta(r-r')>*r^2*dr
    This also calculates the block values of the coordination number
    and the uncertainty according to a 95% confidence interval
    as determined by Student's t-distribution
    """
    nofr= integrate.cumtrapz(gofr*dist**2.*N/L**3.*4*np.pi, dist, initial=0)

    bl_nofr=[] 
    for i in range(1,nblocks+1):
        bl_nofr.append(integrate.cumtrapz(bl_gofr[i-1]*dist**2.*N/L**3.*4*np.pi, dist, initial=0))

    err_nofr = np.std(bl_nofr,axis=0)*t_val # Note: t_val includes the 1/sqrt(nblocks) factor implicitly
    np.savetxt('nofr_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, nofr, err_nofr])
    return nofr, bl_nofr

def calc_dnofr(egofr, dist, bl_egofr, L, N):
    """
    This piece of code calculates the derivative of the coordination number.
    This is achieved via a trapezoid rule integration
    n = -(4*pi*N/V)∫<deltaH delta(r-r')>*r^2*dr
    This also calculates the block values of the coordination number
    and the uncertainty according to a 95% confidence interval
    as determined by Student's t-distribution
    """
    dnofr= integrate.cumtrapz(-egofr*dist**2.*N/L**3.*4*np.pi, dist, initial=0)

    bl_dnofr=[]
    for i in range(1,nblocks+1):
        bl_dnofr.append(integrate.cumtrapz(-bl_egofr[i-1]*dist**2.*N/L**3.*4*np.pi, dist, initial=0))

    err_dnofr = np.std(bl_dnofr,axis=0)*t_val # Note: t_val includes the 1/sqrt(nblocks) factor implicitly
    np.savetxt('dnofr_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, dnofr, err_dnofr])
    return dnofr, bl_dnofr

def pred_nofr(nofr, dnofr, bl_nofr, bl_dnofr, T, Tp):
    """
    This function calculates the first order taylor series prediction based on the nofr and the derivative
    This prediction is saved to a file (i.e. nofr_pred_1_1_235.0_from_298.dat)
    """
    pred = []
    taylorfactor=1/(kb*Tp) - 1/(kb*T) # (beta-beta_o)
    # Makes the prediction for each point in the gofr
    for i in range(len(gofr)):
        pred.append(nofr[i]+dnofr[i]*taylorfactor)
    # Block averaging
    bl_pred = []
    for b in range(nblocks):
        tmp_pred = []
        for i in range(len(bl_nofr[b])):
            tmp_pred.append(bl_nofr[b][i]+bl_dnofr[b][i]*taylorfactor)
        bl_pred.append(tmp_pred)
    err_pred = np.std(bl_pred, axis=0)*t_val
    # Saves the result to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    np.savetxt('nofr_pred_'+str(selec1)+"_"+str(selec2)+"_"+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_pred])
    return pred


def peak_find(dist, gofr, nofr, thresh):
    """
    This is a cute little script for finding the peaks and minima in the RDF
    It prints the value of the maxima, the minima, and the coordination number at that point. 
    The "thresh" parameter determines the threshold by which a peak must stand out compared to
    its neighbors. 
    
    This code writes these values to both the screen, and to a minimum and maximum file.
    """
    # Maxima
    mlocs, dict =  find_peaks(gofr, threshold=thresh, height=1)
    maxvals = []
    dmax = []
    nmax = []
    for i in mlocs:
        print('Maximum value at %4.5f %4.5f with coordination number %4.5f' %(dist[i],gofr[i],nofr[i]))
        nmax.append(nofr[i])
        maxvals.append(gofr[i])
        dmax.append(dist[i])
    minvals = []
    dmin = []
    nmin = []
    for i in range(len(mlocs)-1):
        loc = np.where(gofr == min(gofr[mlocs[i]:mlocs[i+1]]))[0][0]
        minvals.append(gofr[loc])
        dmin.append(dist[loc])
        nmin.append(nofr[loc])
        print('Minimum value at %4.5f %4.5f with coordination number %4.5f' %(dist[loc],gofr[loc],nofr[loc]))
    # Write maximums to file (i.e. gofr_maximums_1_1_298.dat)
    f = open('gofr_maximums_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat','w')
    f.write("#distance(Angstroms) gofrmax nofrmax\n")
    for i in range(len(maxvals)):
        f.write("{0:.5f} {1:.5f} {2:.5f}\n".format(dmax[i],maxvals[i],nmax[i]))
    f.close()
    # Write minimums to file (i.e. gofr_minimums_1_1_298.dat
    g = open('gofr_minimums_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat','w')
    g.write("#distance(Angstroms) gofrmin nofrmin\n")
    for i in range(len(minvals)):
        g.write("{0:.5f} {1:.5f} {2:.5f}\n".format(dmin[i],minvals[i],nmin[i]))
    g.close()
    return maxvals, dmax, minvals, dmin

def calc_PMF(dist, gofr, bl_gofr, T):
    """
    This section of the code calculates the potential of mean force
    The PMF is defined as PMF = -kb*T*log[RDF]
    This is written to a file (i.e. pmf_1_1_298.dat)
    """
    PMF = []
    pdist=[]
    c=0
    for i in range(len(gofr)):
        if i > 1 and gofr[i-1] > 0:
            PMF.append(-kb*T*np.log(gofr[i]))
        else:
            c += 1

    bl_pmf = []
    for b in range(len(bl_gofr)):
        tmp_pmf = []
        for i in range(len(bl_gofr[b])):
            if i > 1 and gofr[i-1] > 0:
                tmp_pmf.append(-kb*T*np.log(bl_gofr[b][i]))
        bl_pmf.append(tmp_pmf)
    err_pmf = np.std(bl_pmf,axis=0)*t_val
    np.savetxt("pmf_"+str(selec1)+"_"+str(selec2)+'_'+str(int(T))+".dat", np.c_[dist[c:], PMF,err_pmf])
    return PMF, c, bl_pmf

def calc_PMF_deriv(PMF, dist, gofr, egofr, bl_gofr, bl_egofr, bl_PMF,T):
    """
    This script calculates the first derivative of the PMF
    The derivative is defined as:
    kb*T*(eRDF/RDF - PMF)
    This is written to a file(i.e. dpmf_1_1_298.dat)
    """
    dPMF = []
    for i in range(len(PMF)):
        dPMF.append(kb*T*(egofr[i]/gofr[i]-PMF[i]))
    bl_dPMF = []
    for b in range(len(bl_PMF)):
        tmp_dPMF = []
        for i in range(len(bl_PMF[b])):
            tmp_dPMF.append(kb*T*(bl_egofr[b][i]/bl_gofr[b][i]-bl_PMF[b][i]))
        bl_dPMF.append(tmp_dPMF)
    err_dPMF = np.std(bl_dPMF,axis=0)*t_val
    np.savetxt("dpmf_"+str(selec1)+"_"+str(selec2)+'_'+str(int(T))+".dat", np.c_[dist, dPMF, err_dPMF])
    return dPMF, bl_dPMF

def pred_PMF(dist, PMF, dPMF,bl_PMF, bl_dPMF, T, Tp):
    """
    This script makes predictions of the PMF at other temperature(s) Tp from the simulation temperature T
    This is done by a standard taylor series expansion out to **one** term
    PMF(beta) = PMF(beta_o) + dPMF/dbeta|beta_o*(beta-beta_o)
    This is written to a file (i.e. pmf_pred_1_1_235.0_from_298.dat)
    """
    pred = []
    taylorfactor=1/(kb*Tp)-1/(kb*T)
    print(taylorfactor)
    for i in range(len(PMF)):
        pred.append(PMF[i]+dPMF[i]*taylorfactor)
    # Block Calculation
    bl_pmfpred = []
    for b in range(len(bl_PMF)):
        tmp_pPMF = []
        for i in range(len(bl_PMF[b])):
            tmp_pPMF.append(bl_PMF[b][i]+bl_dPMF[b][i]*taylorfactor)
        bl_pmfpred.append(tmp_pPMF)
    err_pmfpred = np.std(bl_pmfpred,axis=0)*t_val
    np.savetxt('pmf_pred_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_pmfpred])
    return pred

def calc_dS(PMF, dPMF,bl_PMF,bl_dPMF, dist, T):
    """
    This is a calculation of the entropy from the derivative of the PMF
    dS = 1/(kb*T^2)*(dPMF/dbeta + 2*(kb*T)^2 * log[r])
    Results are written to a file (i.e. dS_1_1.dat)
    """
    dS = []
    for i in range(len(PMF)):
        dS.append(1/(kb*T**2.)*(dPMF[i]+2*(kb*T)**2.*np.log(dist[i])))
    bl_dS = []
    for b in range(len(bl_PMF)):
        tmp_dS = []
        for i in range(len(bl_PMF[b])):
            tmp_dS.append(1/(kb*T**2.)*(bl_dPMF[b][i]+2*(kb*T)**2.*np.log(dist[i])))
        bl_dS.append(tmp_dS)
    err_dS = np.std(bl_dS, axis=0)*t_val
    np.savetxt('dS_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,dS,err_dS])
    return dS



    


import numpy as np
import argparse, math
import matplotlib.pyplot as plt
from scipy import stats,integrate
from scipy.signal import find_peaks

kb=0.0019872041 #kcal/mol


# Command Line Input
parser = argparse.ArgumentParser()
parser.add_argument('-p', default=0, help='[0] No plot [1] Plots in matplotlib')
parser.add_argument('-selec1',default=1, help='Selection1')
parser.add_argument('-selec2', default=1, help='Selection2')
parser.add_argument('-nblocks',default=5, help='Number of blocks')
parser.add_argument('-L', default=21.752, help='Box side length')
parser.add_argument('-N', default=343, help='Number of atoms')
parser.add_argument('-T', default=298.15, help='Temperature')
parser.add_argument('-Tpred', default=[], help='Temperature to predict',action='append')
parser.add_argument('-thresh', help='Threshold for peak finding')
parser.add_argument('-ew', default=0, help='Energy Weighting [0] off [1] on')
parser.add_argument('-graphmovie', default=0, help='Make movie [0] off [1] on')
args = parser.parse_args()
selec1  = str(args.selec1)
selec2  = str(args.selec2)
nblocks = int(args.nblocks)
L       = float(args.L)
N       = int(args.N)
T       = float(args.T)
Tp      = args.Tpred
ew      = int(args.ew)
graphmovie = int(args.graphmovie)
pltflag = int(args.p)


for i in range(len(Tp)):
    Tp[i] = float(Tp[i])

if args.thresh is not None:
    thresh=float(args.thresh)
else:
    thresh=0.001

t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)


# Read GofR
dist, gofr, bl_gofr, err_gofr = read_gofr(t_val,selec1, selec2, nblocks)

# Read dHGofR
if ew == 1:
    egofr, bl_egofr, err_egofr = read_egofr(t_val, selec1, selec2, nblocks)
    pGofR =  []
    if graphmovie == 0:
        for tp in Tp:
            pGofR.append(pred_gofr(gofr, egofr, bl_gofr, bl_egofr, T, float(tp)))
    else:
        for tp in range(220,400):
            pGofR.append(pred_gofr(gofr, egofr, bl_gofr, bl_egofr, T, float(tp)))

# Calculates the Coordination Number
nofr, bl_nofr = calc_nofr(gofr, dist, bl_gofr, L, N)

# Finds the local minima and maxima
maxvals, dmax, minvals, dmin = peak_find(dist, gofr, nofr,thresh)

# Calculates the PMF
PMF, c, bl_PMF = calc_PMF(dist,gofr, bl_gofr, T)


# Calculates the Derivative of the PMF
if ew == 1:
    dnofr, bl_dnofr = calc_dnofr(egofr, dist, bl_egofr, L, N)
    pnofr = []
    for tp in Tp:
        pnofr.append(pred_nofr(nofr,dnofr,bl_nofr, bl_dnofr, T, float(tp)))

    bl_cutgofr = [ bl_gofr[b][c:] for b in range(nblocks) ]
    bl_cutegofr = [ bl_egofr[b][c:] for b in range(nblocks) ]
    dPMF,bl_dPMF=calc_PMF_deriv(PMF, dist[c:], gofr[c:], egofr[c:], bl_cutgofr, bl_cutegofr, bl_PMF, T)

# Makes PMF Prediction
if ew == 1:
    pPMF =  []
    for tp in Tp:
        pPMF.append(pred_PMF(dist[c:], PMF, dPMF,bl_PMF,bl_dPMF, T, tp))

# Calculates Entropy Derivative
if ew == 1:
    dS = calc_dS(PMF, dPMF,bl_PMF,bl_dPMF, dist[c:], T)

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
    



