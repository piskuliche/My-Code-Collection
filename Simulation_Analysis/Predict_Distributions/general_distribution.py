#!/usr/bin/env python
import numpy as np
import pickle
import os,time, argparse
from scipy import stats
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import pdist, squareform
from asphere import wrap_box, polyhedron, compute_vc, asphericity

def pbc(r1,r2):
    """ 
    This calculates the periodic boundary conditions minimum image distance between two points.
    """
    dr = np.subtract(r1,r2)
    #dist = np.linalg.norm(dr - np.round(dr))
    dist = np.sqrt(np.sum(np.subtract(dr,np.round(dr))**2.))
    return dist

def pbc2(r1,r2,otmp):
    """
    This calculates the PBC minimum image distance between two numpy arrays, of indeterminant length
    and then returns the minimum of all the minimum image distances
    """
    dr = np.subtract(r1,r2)
    dist = np.sqrt(np.sum(np.subtract(dr,np.round(dr))**2.,axis=1))
    #dist = np.linalg.norm(np.subtract(dr,np.round(dr)),axis=1)
    minval = dist.min()
    atom = otmp[np.where(dist==minval)]
    return minval,atom


def find_closest(dr_arr, drop):
    """
    This finds the minimum values
    """
    m,n=dr_arr.shape
    mask=np.ones((m,n),dtype=bool)
    mask[range(m),drop]=False
    finarray=dr_arr[mask].reshape(m,n-1)
    minval=finarray.min(axis=1)
    minh,mino=np.where(dr_arr==finarray.min(axis=1,keepdims=True))
    return minval, mino
    
    

def distance_wrap(r):
    """
    This calculates the distances with periodic boundaries
    """
    dr=np.sqrt(squareform(pdist(r[:,0])-np.round(pdist(r[:,0])))**2. + squareform(pdist(r[:,1])-np.round(pdist(r[:,1])))**2. + squareform(pdist(r[:,2])-np.round(pdist(r[:,2])))**2.)
    return dr

def calc_hbonds(frame):
    """
    This calculates the oo, oh, and angle distributions for hydrogen bond exchanges.
    """
    roo,roh,cosang=[],[],[]
    natoms = len(frame["type"])
   
    # Masks dr to be based on hatoms and all others
    hdr=np.array(frame["dr"][frame["h"]])
    hodr=np.array(hdr[:,frame["o"]])
    # Masks dr ot be based on oatoms 
    odr=np.array(frame["dr"][frame["o"]])

    # This is an array that is 686 long and tells what  each hatom is closest to
    oatms=(np.array(frame["co"])/3)[frame["h"]].astype(int)

    # Calculates distances and angle.
    sides_roh,closest=find_closest(hodr,oatms)
    sides_rho=hodr.min(axis=1)
    oarr=odr[:,np.array(frame["o"])]
    sides_roo=oarr[oatms,closest]
    roh = np.multiply(sides_roh,frame["L"])
    rho = np.multiply(sides_rho,frame["L"])
    roo = np.multiply(sides_roo,frame["L"])
    #Law of cosines to get the hbond angle
    #cos_theta=(side_rho**2.+side_roo**2.-side_roh**2.)/(2*side_rho*side_roo)
    cosang = np.divide(np.power(rho,2)+np.power(roo,2)-np.power(roh,2),np.multiply(2,np.multiply(rho,roo)))
    # Does histogramming of distances and angles
    histoh,bins=np.histogram(roh,bins=50,range=(1.0,4.0),density=False)
    histoh = histoh/len(roh)
    histoo,bins=np.histogram(roo,bins=50,range=(1.0,4.0),density=False)
    histoo = histoo/len(roo)
    histang, bins=np.histogram(np.arccos(cosang), bins=50, range=(0.0,1.0), density=False)
    histang = histang/len(cosang)
    return np.asarray(histoo),np.asarray(histoh),np.asarray(histang)


def do_analysis(params,frame):
    """ 
    This allows for modular input and output, and provides timings.
    """
    roo,roh,ctheta,asp=np.zeros(50),np.zeros(50),np.zeros(50),np.zeros(50)
    roo,roh,ctheta=calc_hbonds(frame)
    if params["aspon"]==1: asp=asphericity(frame)
    return roo, roh,ctheta, asp

def write_data(params,finaldata,energy):
    """
    This writes the data to an outut file
    """
    # Does all the zeroing
    rbins=np.linspace(1.03,4.03,50)
    cbins=np.linspace(0.1,1.05,50)
    for key in finaldata:
        if 'err' not in key and "bl" not in key:
            if 'C' not in key:
                np.savetxt(params["pre"]+key+'_dist.dat', np.c_[rbins,finaldata[key],finaldata[key+'err']])
            else:
                np.savetxt(params["pre"]+key+'_dist.dat', np.c_[cbins,finaldata[key],finaldata[key+'err']])
    return

def post_analysis(params,postdata,energy):
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

        

def read_traj(params):
    """
    This reads in an xyz trajectory file, and then calculates the derivatives, etc
    """
    natoms=0
    # Opens the file to read number of atoms
    with open(params["filename"]) as f:
        line=f.readline().strip()
        natoms=int(line)
    # Read in the distances
    print("Reading Box Length")
    Lval=np.genfromtxt('L.dat',usecols=0,max_rows=params["stop"])
    # Read in the energy, and calculate the fluctuation of energy
    print("Reading in Energies")
    e,lj,ke,vol = [], [], [], []
    if os.path.isfile("e_init.out"):   e=np.genfromtxt('e_init.out',usecols=0,max_rows=params["stop"])[:params["stop"]]
    if os.path.isfile("lj_init.out"):  lj=np.genfromtxt('lj_init.out',usecols=0,max_rows=params["stop"])[:params["stop"]]
    if os.path.isfile("vol_init.out"): vol=np.genfromtxt('vol_init.out',usecols=0,max_rows=params["stop"])[:params["stop"]]
    if os.path.isfile("ke_init.out"):  ke=np.genfromtxt('ke_init.out',usecols=0,max_rows=params["stop"])[:params["stop"]]
    energy = { "e":e, "lj":lj, "ke":ke, "vol":vol }
    if len(lj) > 0 and len(ke) > 0 and len(lj) > 0:
        elec = np.subtract(e,lj)
        elec = np.subtract(elec,ke)
        energy["elec"]=elec

    # Here we define hs and os
    print("Defining Essential Parameters")
    h,o,co,mol,atype=[],[],[],[],[]
    if params["ovtype"]=="ohh":
        count,t=0,0
        for i in range(natoms):
            if i%3 == 0:
                count+=1
                o.append(i)
                co.append(i)
                atype.append(1)
                t=i
            else:
                h.append(i)
                atype.append(2)
                co.append(t)
            mol.append(count)
    data={ "Roo":[], "Roh":[], "Ctheta":[], "Asphere":[]}
    # Starts reading the file again, this time for real
    print("Opening Trajectory File")
    with open(params["filename"]) as f:
        lperframe=natoms+2
        frame=[]
        framecount=0
        # Each iteration is a loop over frames.
        start = time.time()
        while True:
            # Creates dictionary
            frame={ "type":atype,"co":co, "r":[],"ra":[],"mol":mol, "h":h, "o":o,"L":Lval[framecount]}
            # Skips the initial two lines of the xyz file
            line=f.readline()
            line=f.readline()
            if not line:
                break
            # Reads in the frame into the dictionary, "frame"
            for l in range(lperframe-2):
                line=f.readline()
                vals=line.strip().split()
                frame["ra"].append(np.array((float(vals[1])/frame["L"],float(vals[2])/frame["L"],float(vals[3])/frame["L"])))
                frame["r"].append([[float(vals[1])/frame["L"]],[float(vals[2])/frame["L"]],[float(vals[3])/frame["L"]]])
                # increments mol number every 3 atoms
            # Do analysis
            frame["r"]=np.array(frame["r"])
            frame["dr"]=distance_wrap(frame["r"])
            roo,roh,ctheta,asp=do_analysis(params,frame)
            # The following section increments storage vectors
            data["Roo"].append(roo)
            data["Roh"].append(roh)
            data["Ctheta"].append(ctheta*180.0/np.pi)
            data["Asphere"].append(asp)
            # End storage Vector section
            framecount+=1
            end = time.time()
            if (framecount % params["stop"] == 0): break
            if (framecount % 10 == 0): print("frame: %s \ntime_per_frame: %s seconds\ntotal_time: %s seconds" % (framecount,(end-start)/framecount, end-start))
            # Writes a restart file
            if (framecount % 1000 == 0):
                g = open("restart_distribution.pkl","wb")
                pickle.dump(data,g)
                g.close()
        manipulate_data(params,data,energy)

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


    

parser = argparse.ArgumentParser()
parser.add_argument('-f', default="traj.xyz", help='Trajectory file name')
parser.add_argument('-nblocks', default=5, help='Number of blcoks')
parser.add_argument('-nconfigs', default=1000, help='Total number of configurations')
parser.add_argument('-oatom', default=1, help='Integer type representing oxygen')
parser.add_argument('-hatom', default=2, help='Integer type representing hydrogen')
parser.add_argument('-order', default="ohh", help='Order of trajectory file, default is ohh')
parser.add_argument('-T', default=298.15, help='Temperature of the simulation (K)')
parser.add_argument('-P', default=1.0, help='Pressure of the simulation (bar)')
parser.add_argument('-prepend', default="run_", help='Prepend the output files with information')
parser.add_argument('-asphere', default=0, help='Boolean value, 1 to calculate asphericity, 0 to not')

args = parser.parse_args()

kb=0.0019872041

inputparams={"filename":str(args.f), "stop":int(args.nconfigs), "htype":int(args.hatom), "otype":int(args.oatom), "ovtype":str(args.order), "nblocks":int(args.nblocks), "T":float(args.T), "P":float(args.P),"pre":str(args.prepend), "aspon":int(args.asphere)}
print("Welcome to the Distribution Predictor!")
if(inputparams["aspon"]==1): print("Note: Asphericity calculation is on, calcualtion will be much slower.")
read_traj(inputparams)
