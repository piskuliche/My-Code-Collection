#!/usr/bin/env python
"""
This code calculates the shear viscosity through the Green-Kubo relation.

eta=V/kbT int_0^t {<P_ab(0)P_ab(t)>dt}

Copyright 2020 Ezekiel Piskulich, University of Kansas

"""

import numpy as np
import argparse
from scipy import stats

consts={"kb":1.370649E-23,"atmtobar":1.01325,"UnitstocentiPoise":1E-32,"UnitstocPperps":1E-29}


def read_log(inputparams):
    """
    Extracts the data from the lammps log file and puts it into a data structure.
    """
    extract=["Step","Volume","Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"]
    data={}
    for item in extract:
        data[item]=[]
    dataloc={}
    flag=0
    with open(inputparams["logfile"],'r') as f:
        lines=[]
        for line in f:
            if "Loop" in line:
                flag=0
            if flag==1:
                for item in extract:
                    data[item].append(float(line.split()[dataloc[item]]))
            if "Step" in line:
                for item in extract:
                    if item in line:
                        dataloc[item]=int(line.split().index(item))
                        flag=1
                    else: 
                        print("ERROR: %s not written to logfile" % item)
                        exit()
    for item in extract:
        data[item]=data[item][inputparams["skip"]:]
    data["nopoints"]=len(data["Pxx"])
    data["dumpfreq"]=float(data["Step"][1]-data["Step"][0])
    return data

def av_press(logdata,to,t):
    """
    Calculates the average pressure contributions.
    """
    Pxy=logdata["Pxy"][to]*logdata["Pxy"][t]
    Pyz=logdata["Pyz"][to]*logdata["Pyz"][t]
    Pxz=logdata["Pxz"][to]*logdata["Pxz"][t]
    Pxxyy=((logdata["Pxx"][to]-logdata["Pyy"][to])/2.)*((logdata["Pxx"][t]-logdata["Pyy"][t])/2.)
    Pyyzz=((logdata["Pyy"][to]-logdata["Pzz"][to])/2.)*((logdata["Pyy"][t]-logdata["Pzz"][t])/2.)
    return (Pxy + Pyz + Pxz + Pxxyy + Pyyzz)/5.


def calc_eta(params,logdata,to):
    """
    Calculates the shear viscosity over time
    """
    eta,avP,tvals = [],[],[]
    prefactor=params["V"]/(consts["kb"]*params["T"])
    for t in range(to,to+params["length"]):
        tvals.append(float(t-to)*logdata["dumpfreq"]*params["dt"])
        avP.append(prefactor*av_press(logdata,to,t))
        eta.append(np.trapz(avP,tvals))
    return eta,avP

def calculate_visc(params):
    """
    The actual code that does the calculation.
    """
    print("Starting calculation")
    logdata=read_log(inputparams)
    print("Log data has been read")
    params["V"]=logdata["Volume"][0]
    total_to=int((logdata["nopoints"]-params["length"])/params["tosep"])
    final_to=int(logdata["nopoints"]-params["length"])
    print("There will be %d time origins and the last will be %d" % (total_to, final_to))
    finaldata={"eta":[],"avP":[]}
    for to in range(0,final_to,int(params["tosep"])):
        if(to%1000 == 0): print(to)
        step_eta,step_avP=calc_eta(params,logdata,to)
        finaldata["eta"].append(step_eta)
        finaldata["avP"].append(step_avP)
    times=np.linspace(0,params["length"]*logdata["dumpfreq"]*params["dt"],params["length"])
    # Write Average Files
    av_eta=np.average(finaldata["eta"],axis=0)*consts["atmtobar"]**2.*consts["UnitstocentiPoise"]
    av_betaV_avP=np.average(finaldata["avP"],axis=0)*consts["atmtobar"]**2.*consts["UnitstocPperps"]
    nperb=int(total_to/params["nblocks"])
    bl_av_eta,bl_av_betaV_avP=[],[]
    for b in range(params["nblocks"]):
        bstart=b*nperb
        bend=(b+1)*nperb
        bl_av_eta.append(np.average(finaldata["eta"][bstart:bend],axis=0)*consts["atmtobar"]**2.*consts["UnitstocentiPoise"])
        bl_av_betaV_avP.append(np.average(finaldata["avP"][bstart:bend],axis=0)*consts["atmtobar"]**2.*consts["UnitstocPperps"])
        np.savetxt("bl_"+str(b)+"_shear.dat", np.c_[times,bl_av_eta[b],bl_av_betaV_avP[b]])
    std_eta=np.std(bl_av_eta,axis=0)*params["tval"]
    std_betaV_avP=np.std(bl_av_betaV_avP,axis=0)*params["tval"]
    np.savetxt("shear.dat", np.c_[times,av_eta,std_eta,av_betaV_avP,std_betaV_avP])
    return


parser = argparse.ArgumentParser()
parser.add_argument('-f', default="log.lammps",type=str, help='Log file to parse')
parser.add_argument('-skip', default=0, type=int, help='Number of lines to skip')
parser.add_argument('-nb', default=5, type=int, help='Number of blocks')
parser.add_argument('-to', default=10, type=int, help='Number of steps between origins')
parser.add_argument('-len', default=500, type=int, help='Number of steps to calculate each corr_func')
parser.add_argument('-T', default=298.15, type=float, help='Temperature of simulation')
parser.add_argument('-timestep', default=1.0, type=float, help='Timestep in fs')
args= parser.parse_args()

inputparams={"logfile":args.f, "skip":args.skip, "nblocks":args.nb,"length":args.len, "tosep":args.to, "T":args.T, "dt":args.timestep}
inputparams["tval"]=stats.t.ppf(0.975,inputparams["nblocks"]-1)/np.sqrt(inputparams["nblocks"])

calculate_visc(inputparams)
