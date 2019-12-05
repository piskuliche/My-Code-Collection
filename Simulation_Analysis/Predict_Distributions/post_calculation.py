import numpy as np
import pickle
import os,time, argparse
from scipy import stats
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import pdist, squareform

kb=0.0019872041


def post_analysis(params, postdata, energy):
    """
    This calls all the output routines and does the calculations
    """
    t_val=stats.t.ppf(0.975,params["nblocks"]-1)/np.sqrt(params["nblocks"])
    post_out={}
    for key in postdata:
        if 'err' not in key and "bl" not in key and "e" not in key[:2] and "lj" not in key and "vol" not in key:
            print(key)
            post_out["AG"+key],post_out["UH"+key],post_out["S"+key]=calc_thermodynamic_potential(params, postdata[key], postdata["e"+key])
            tmpag, tmpuh, tmps = [], [], []
            for b in range(params["nblocks"]):
                ag,uh,s=calc_thermodynamic_potential(params, postdata["bl_"+key][b], postdata["bl_e"+key][b])
                tmpag.append(ag)
                tmpuh.append(uh)
                tmps.append(s)
            post_out["AG"+key+"err"] = np.std(tmpag,axis=0)*t_val
            post_out["UH"+key+"err"] = np.std(tmpuh,axis=0)*t_val
            post_out["S"+key+"err"] = np.std(tmps,axis=0)*t_val
            for ekey in energy:
                post_out["UH"+ekey+key]=calc_H_or_U(params, postdata[key], postdata[ekey+key])
                tmpuh=[]
                for b in range(params["nblocks"]):
                    tmpuh.append(calc_H_or_U(params, postdata["bl_"+key][b], postdata["bl_"+ekey+key][b]))
                post_out["UH"+ekey+key+"err"] = np.std(tmpuh,axis=0)*t_val
    return post_out

def write_data(params,finaldata,energy):
    """
    This writes the data to an outut file
    """
    # Does all the zeroing
    rmin,rmax=1.0,4.0
    cmin,cmax=0.0,1.0
    rbin,cbin=200.,200.
    rdiff=(rmax-rmin)/(2*rbin)
    cdiff=(rmax-rmin)/(2*cbin)
    rbins=np.linspace(rmin+rdiff,rmax-rdiff,rbin)
    cbins=np.linspace(cmin+cdiff,cmax-cdiff,cbin)
    for key in finaldata:
        if 'err' not in key and "bl" not in key:
            if 'C' not in key:
                np.savetxt(params["pre"]+key+'_dist.dat', np.c_[rbins,finaldata[key],finaldata[key+'err']])
            else:
                np.savetxt(params["pre"]+key+'_dist.dat', np.c_[cbins,finaldata[key],finaldata[key+'err']])
    return

def manipulate_data(params, data, energy):
    """
    This code takes in the energy files and averages
    """
    t_val=stats.t.ppf(0.975,params["nblocks"]-1)/np.sqrt(params["nblocks"])
    # Checks if energies have data in them
    for key in energy:
        if len(energy[key]) == 0: del energy[key]
    # Create Output Data Structure
    OutputData = {}
    nperb = int(params["stop"]/params["nblocks"])
    # Total Calculation
    for key in data:
        OutputData[key] = np.average(data[key],axis=0)
        tmp = []
        for b in range(params["nblocks"]):
            bstart, bend = b*nperb, (b+1)*nperb
            tmp.append(np.average(data[key][bstart:bend],axis=0))
        OutputData[key+"err"]=np.std(tmp,axis=0)*t_val
        OutputData["bl_"+key]=tmp
        for ekey in energy:
            OutputData[ekey+key]=np.average(np.array(data[key])*np.array(energy[ekey])[:,None]-OutputData[key]*np.average(energy[ekey]),axis=0)
            tmp = []
            for b in range(params["nblocks"]):
                bstart, bend = b*nperb, (b+1)*nperb
                tmp.append(np.average(np.array(data[key][bstart:bend])*np.array(energy[ekey][bstart:bend])[:,None],axis=0)-OutputData[key]*np.average(energy[ekey][bstart:bend]))
            OutputData[ekey+key+"err"]=np.std(tmp,axis=0)*t_val
            OutputData["bl_"+ekey+key]=tmp
    write_data(params,OutputData,energy)
    PostData=post_analysis(params,OutputData,energy)
    write_data(params,PostData,energy)
    return

def calc_thermodynamic_potential(params,probdist,eprobdist):
    """
    This calculates the free energy from a probability distribution.
    A(NVT)=-kb*T*np.log(P)
    G(NPT)=-kb*T*np.log(P)
    """
    AG = -kb*params["T"]*np.log(probdist,out=np.zeros(len(probdist)), where=probdist!=0)
    UH = calc_H_or_U(params,probdist,eprobdist)
    S  = calc_S(params,AG, UH)
    return AG, UH, S

def calc_H_or_U(params, probdist, eprobdist):
    """
    Calculates the enthalpy or internal energy
    U(NVT)=Ph/P
    H(NPT)=Ph/P
    """
    UH=np.divide(eprobdist,probdist,out=np.zeros(len(eprobdist)),where=probdist!=0)
    return UH

def calc_S(params, AG,UH):
    """
    Calculates the entropy
    S = (UH-AG)/T
    """
    S = (UH-AG)/params["T"]
    return S
