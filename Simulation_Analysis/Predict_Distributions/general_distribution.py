#!/usr/bin/env python
import numpy as np
import time, argparse
from scipy import stats


def pbc(r1,r2):
    """ 
    This calcualtes the periodic boundary conditions minimum image distance between two points.
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

def calc_hbonds(frame):
    """
    This calculates the oo, oh, and angle distributions for hydrogen bond exchanges.
    """
    roo,roh,cosang=[],[],[]
    natoms = len(frame["type"])
    closest=[0,0]
    for atom1 in frame["h"]:
        #oatm=frame["co"][atom1]
        oatm=int(frame["co"][atom1]/3)
        otmp = np.delete(frame["o"],oatm)
        closest=pbc2(frame["r"][atom1],np.array(frame["r"])[otmp],otmp)
        # Scans other oxygens for closest oxygen to current selection
        # Saves oxygen oxygen distance for histogramming
        side_roh=closest[0]*frame["L"]
        roh.append(side_roh)
        # Identifies hydrogen in hbond
        partner=closest[1][0]
        side_roo = pbc(frame["r"][frame["o"][oatm]],frame["r"][partner])*frame["L"]
        roo.append(side_roo)
        side_rho=np.hypot(side_roo,side_roh)
        # Calculate the cosine of the angle
        # cos(theta)=(rho^2+roo^2-roh^2)/(2*rho*roo)
        #print(side_rho,side_roo,side_roh)
        cos_theta=(side_rho**2.+side_roo**2.-side_roh**2.)/(2*side_rho*side_roo)
        cosang.append(cos_theta)
    #roo,roh,cosang = zip(*Parallel(n_jobs=4)(delayed(process_it)(i,frame) for i in frame["h"]))
    # Does histogramming of distances and angles
    print(np.array(roh).max())
    histoh,bins=np.histogram(roh,bins=100,range=(0.0,5.0),density=False)
    histoh = histoh/len(roh)
    histoo,bins=np.histogram(roo,bins=100,range=(0.0,5.0),density=False)
    histoo = histoo/len(roo)
    histang, bins=np.histogram(cosang, bins=100, range=(0.0,1.0), density=False)
    histang = histang/len(cosang)
    return np.asarray(histoo),np.asarray(histoh),np.asarray(histang)


def do_analysis(frame):
    """ 
    This allows for modular input and output, and provides timings.
    """
    start=time.time()
    roo,roh,ctheta=calc_hbonds(frame)
    end=time.time()
    print(end-start)
    return roo, roh,ctheta

def write_data(params,finaldata,energy):
    """
    This writes the data to an outut file
    """
    #np.savetxt('r_oo_vs_r_oh.dat',np.c_[rbins,Roo/float(framecount),Roh/float(framecount),(dRoo)/float(framecount),(dRoh)/float(framecount)])
    #np.savetxt('c_theta.dat', np.c_[cbins, Ctheta/float(framecount), dCtheta/float(framecount)])
    for key in finaldata:
        np.savetxt(key+'_dist.dat', np.c_[rbins,cbins,finaldata[key]],finaldata[key+'err'])
        for ekey in energy:
            np.savetxt(ekey+key+'_dist.dat', np.c_[rbins,cbins,finaldata[ekey+key]],finaldata[ekey+key+'err']])
    return

def manipulate_data(params, data, energy):
    """
    This code takes in the energy files and averages
    """
    t_val=stats.t.ppf(0.975,params["nblocks"]-1)/np.sqrt(params["nblocks"])
    # Checks if energies have data in them
    for key in energy:
        if len(energy[key]) = 0: del energy[key]
    # Create Output Data Structure
    OutputData = {}
    nperb = int(params["nconfigs"]/params["nblocks"])
    # Total Calculation
    for key in data:
        OutputData[key] = np.average(data[key],axis=0)
        tmp = []
        for b in range(params["nblocks"]):
            bstart, bend = b*nperb, (b+1)*nperb
            tmp.append(np.average(data[key][bstart:bend],axis=0))
        OutputData[key+"err"]=np.std(tmp)*t_val
        for ekey in energy:
            OutputData[ekey+key]=np.average(np.array(data[key])*np.array(energy[ekey])[:,None]-OutputData[key]*np.average(energy[ekey]).axis=0)
            tmp = []
            for b in range(params["nblocks"]):
                bstart, bend = b*nperb, (b+1)*nperb
                tmp.append(np.average(np.array(data[key][bstart:bend])*np.array(energy[ekey])[:,None]-OutputData[key]*np.average(energy[ekey][bstart:bend]),axis=0))
            OutputData[ekey+key+"err"]=np.std(tmp)*t_val
    write_data(params,finaldata,energy)
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
    Lval=np.genfromtxt('L.dat',usecols=0)
    # Read in the energy, and calculate the fluctuation of energy
    e,lj,ke,vol = [], [], [], []
    if os.isfile("e_init.out"):   e=np.genfromtxt('e_init.out',usecols=0)[:stop]
    if os.isfile("lj_init.out"):  lj=np.genfromtxt('lj_init.out',usecols=0)[:stop]
    if os.isfile("vol_init.out"): vol=np.genfromtxt('vol_init.out',usecols=0)[:stop]
    if os.isfile("ke_init.out"):  ke=np.genfromtxt('ke_init.out',usecols=0)[:stop]
    energy = { "e":e, "lj":lj, "ke":ke, "vol":vol }
    # Does all the zeroing
    rbins=np.linspace(0.025,5.025,100)
    cbins=np.linspace(0.005,1.005,100)
    # Here we define hs and os
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
    data={ "Roo"=[], "Roh"=[], "Ctheta"=[] }
    # Starts reading the file again, this time for real
    with open(params["filename"]) as f:
        lperframe=natoms+2
        frame=[]
        framecount=0
        # Each iteration is a loop over frames.
        while True:
            # Creates dictionary
            frame={ "type":atype,"co":co, "r":[],"mol":mol, "h":h, "o":o,"L":Lval[framecount]}
            # Skips the initial two lines of the xyz file
            line=f.readline()
            line=f.readline()
            if not line:
                break
            # Reads in the frame into the dictionary, "frame"
            for l in range(lperframe-2):
                line=f.readline()
                vals=line.strip().split()
                frame["r"].append(np.array((float(vals[1])/frame["L"],float(vals[2])/frame["L"],float(vals[3])/frame["L"])))
                # increments mol number every 3 atoms
            # Do analysis
            roo,roh,ctheta=do_analysis(frame)
            # The following section increments storage vectors
            data["Roo"].append(roo)
            data["Roh"].append(roh)
            data["Ctheta"].append(ctheta)
            # End storage Vector section
            framecount+=1
            if (framecount % params["nconfigs"] == 0): break
            if (framecount % 100 == 0): print(framecount)
        manipulate_data(params,data,energy)
    return

parser = argparse.ArgumentParser()
parser.add_argument('-f', default="traj.xyz", help='Trajectory file name')
parser.add_argument('-nblocks', default=5, help='Number of blcoks')
parser.add_argument('-nconfigs', default=1000, help='Total number of configurations')
parser.add_argument('-o', default=1, help='Integer type representing oxygen')
parser.add_argument('-h', default=2, help='Integer type representing hydrogen')
parser.add_argument('-order', default="ohh", help='Order of trajectory file, default is ohh')

args = parser.parse_args()


inputparams={"filename":str(args.f), "stop"=int(args.nconfigs), "htype":int(args.h), "otype":int(args.o), "ovtype":str(args.order)}
read_traj(inputparams)
