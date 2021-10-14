#!/usr/bin/env python
import numpy as np
from scipy import stats
import sys

nbins = int(sys.argv[1])
natms = int(sys.argv[2])
nsubs = int(sys.argv[3])
nblocks = int(sys.argv[4])

t_val = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

def hist_data(filename,nbins,counts):
    x, hist, norm = [],[],[]
    for subdir in range(1,nsubs):
        tmpx, tmphist, tmpnorm = np.genfromtxt(str(subdir)+"/"+filename, usecols=(0,1,2),dtype=float, unpack=True)
        x=np.append(x, tmpx)
        hist=np.append(hist,tmphist)
        norm=np.append(norm, tmpnorm)

    # Reshape hist by # bin axis
    norm = np.reshape(norm, (-1, nbins))
    hist = np.reshape(hist, (-1, nbins))
    singleX = np.reshape(x, (-1, nbins))[0]
    if counts == 0:
        counts = np.sum(hist,axis=1)[0]
    #totalhist=np.sum(np.divide(hist.T,counts),axis=1)
    totalhist = np.sum(np.divide(hist.T, norm.T), axis=1)
    # Do the blocking
    bl_hist_raw=np.array_split(hist,nblocks)
    norm_raw=np.array_split(norm,nblocks)
    bl_hist=[]
    for bl in range(nblocks):
        bl_hist.append(np.sum(np.divide(np.array(bl_hist_raw[bl]).T,np.array(norm_raw[bl]).T),axis=1))
    error=np.std(bl_hist,axis=0)*t_val
    np.savetxt("final_"+filename, np.c_[singleX,totalhist/counts,error/counts])
    return counts

def thick_data(filename):
    hdr, ldr, rvc = [],[],[]
    for subdir in range(1,nsubs):
        tmphdr,tmpldr,tmprvc = np.genfromtxt(str(subdir)+"/"+filename, usecols=(1,2,3), unpack=True)
        hdr=np.append(hdr,tmphdr)
        ldr=np.append(ldr,tmpldr)
        rvc=np.append(rvc,tmprvc)
    # Calculate Block Values
    bl_hdr = np.average(np.array_split(hdr,nblocks),axis=1)
    bl_ldr = np.average(np.array_split(ldr,nblocks),axis=1)
    bl_rvc = np.average(np.array_split(rvc,nblocks),axis=1)
    # Calculate the Error
    error_hdr = np.std(bl_hdr)*t_val
    error_ldr = np.std(bl_ldr)*t_val
    error_rvc = np.std(bl_rvc)*t_val
    # Calculate the Surface Area
    SA = 4 * np.pi * hdr**2
    bl_SA = 4*np.pi*np.power(bl_hdr,2)
    error_SA = np.std(bl_SA)*t_val
    print("Outer Shell Distance from COM")
    print("%10.5f +/- %10.5f" % (np.average(hdr),error_hdr))
    print("Inner Shell Distance from COM")
    print("%10.5f +/- %10.5f" % (np.average(ldr),error_ldr))
    print("Vesicle Radius")
    print("%10.5f +/- %10.5f" % (np.average(rvc), error_rvc))
    print("Surface Area")
    print("%10.5f +/- %10.5f" % (np.average(SA), error_SA))
    f=open("finalvalues_"+filename,'w')
    f.write("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %(np.average(hdr),error_hdr, np.average(ldr),error_ldr,np.average(rvc),error_rvc))
    f.close()

print("Starting")
cnts=hist_data("dr_out_hist_raw.dat", nbins, 0)
thick_data("drh_out.dat")
hist_data("hdr_out_hist_raw.dat", nbins, cnts)
hist_data("ldr_out_hist_raw.dat", nbins, cnts)
for atm in range(1,natms+1):
    print("Atom:",atm)
    hist_data("atomdr_"+str(atm)+"_hist_raw.dat", nbins, cnts)
