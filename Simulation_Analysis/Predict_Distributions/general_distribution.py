#!/usr/bin/env python
import numpy as np
import os,time, argparse
from scipy import stats
from scipy.spatial import ConvexHull, Voronoi




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

def wrap_box(r1, r2):
    """
    This wraps the coordinates (in box units) back into the box.
    """
    dr = np.subtract(r1,r2)
    new_coords = np.subtract(dr,np.round(dr))
    return new_coords

def polyhedron(coords, atom):
    """
    This finds the polyhedra for the center molecule
    """
    vor = Voronoi(coords)
    points = [vor.vertices[x] for x in vor.regions[vor.point_region[int(atom/3)]] if x != -1]
    return points

def compute_vc(points):
    """
    Computes the voronoi cell
    """
    S = ConvexHull(points).area
    V = ConvexHull(points).volume
    # Voronoi cell
    eta = S**3./(36.*np.pi*V**2.)
    return eta


def asphericity(frame):
    """
    Note: this implementation is a hacked together implementation of the
    asphericity calculation as was included in the Iorder package,
    which can be found at https://github.com/ipudu/order/blob/master/order/avc.py

    """
    e=[]
    for atom1 in frame["o"]:
        c = frame["r"][atom1]
        cs = np.asarray(frame["r"])[np.asarray(frame["o"])]
        nc = wrap_box(c, cs)
        points = polyhedron(nc, atom1)
        e.append(compute_vc(points))
    histasp,bins =np.histogram(e,bins=50,range=(1.0,4.0),density=False)
    histasp = histasp/len(e)
    return np.asarray(histasp)


def calc_hbonds(frame):
    """
    This calculates the oo, oh, and angle distributions for hydrogen bond exchanges.
    """
    roo,roh,cosang=[],[],[]
    natoms = len(frame["type"])
    closest=[0,0]
    for atom1 in frame["h"]:
        # Bookkeeping
        # atom1 is the donor hydrogen atom
        # oatm is its donor oxygen atom
        # partner is the oxygen acceptor
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
        #print(frame["L"])
        #print("hd",frame["r"][atom1]*frame["L"])
        #print("od",frame["r"][frame["o"][oatm]]*frame["L"])
        #print("oa",frame["r"][partner]*frame["L"])
        #side_rho=np.sqrt(np.abs(side_roo**2.-side_roh**2.))
        side_rho=pbc(frame["r"][atom1],frame["r"][frame["o"][oatm]])*frame["L"]
        #print(side_roh, side_roo, side_rho)
        #side_rho=0.9572
        # Calculate the cosine of the angle
        # cos(theta)=(rho^2+roo^2-roh^2)/(2*rho*roo)
        cos_theta=(side_rho**2.+side_roo**2.-side_roh**2.)/(2*side_rho*side_roo)
        cosang.append(cos_theta)
    #roo,roh,cosang = zip(*Parallel(n_jobs=4)(delayed(process_it)(i,frame) for i in frame["h"]))
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
    #np.savetxt('r_oo_vs_r_oh.dat',np.c_[rbins,Roo/float(framecount),Roh/float(framecount),(dRoo)/float(framecount),(dRoh)/float(framecount)])
    #np.savetxt('c_theta.dat', np.c_[cbins, Ctheta/float(framecount), dCtheta/float(framecount)])
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
    Lval=np.genfromtxt('L.dat',usecols=0)
    # Read in the energy, and calculate the fluctuation of energy
    e,lj,ke,vol = [], [], [], []
    if os.path.isfile("e_init.out"):   e=np.genfromtxt('e_init.out',usecols=0)[:params["stop"]]
    if os.path.isfile("lj_init.out"):  lj=np.genfromtxt('lj_init.out',usecols=0)[:params["stop"]]
    if os.path.isfile("vol_init.out"): vol=np.genfromtxt('vol_init.out',usecols=0)[:params["stop"]]
    if os.path.isfile("ke_init.out"):  ke=np.genfromtxt('ke_init.out',usecols=0)[:params["stop"]]
    energy = { "e":e, "lj":lj, "ke":ke, "vol":vol }
    if len(lj) > 0 and len(ke) > 0 and len(lj) > 0:
        elec = np.subtract(e,lj)
        elec = np.subtract(elec,ke)
        energy["elec"]=elec

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
    data={ "Roo":[], "Roh":[], "Ctheta":[], "Asphere":[]}
    # Starts reading the file again, this time for real
    with open(params["filename"]) as f:
        lperframe=natoms+2
        frame=[]
        framecount=0
        # Each iteration is a loop over frames.
        start = time.time()
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
read_traj(inputparams)
