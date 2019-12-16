#!/usr/bin/env python
"""
This software was developed by Ezekiel Piskulich at the University of Kansas in November 2019.

Feel free to modify, extend, or develop this code further.

"""




import numpy as np
import pickle
import os,time, argparse
from scipy import stats
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import pdist, squareform
from asphere import wrap_box, polyhedron, compute_vc, asphericity
from rdf import radial_distribution, remove_likes, init_rdf 
from hbond_calc import find_closest, distance_wrap, calc_hbonds
from post_calculation import calc_S, calc_H_or_U, calc_thermodynamic_potential
from post_calculation import manipulate_data, write_data, post_analysis
from restart import call_calculation


def do_analysis(params,frame):
    """ 
    This allows for modular input and output, and provides timings.
    """
    roo,roh,ctheta,asp=np.zeros(200),np.zeros(200),np.zeros(200),np.zeros(200)
    goo, goh, ghh = np.zeros(params['rdfbins']),np.zeros(params['rdfbins']),np.zeros(params['rdfbins'])
    roo,roh,ctheta=calc_hbonds(frame)
    if params["rdfon"]==1: goo,goh,ghh = radial_distribution(params,frame)
    if params["aspon"]==1: asp=asphericity(frame)
    return roo, roh, ctheta, asp, goo, goh, ghh

def farm_frames(params):
    print("Preparing to farm out frames")
    natoms=0
    # Opens the file to read number of atoms
    with open(params["filename"]) as f:
        line=f.readline().strip()
        natoms=int(line)
    print("There are %d atoms per frame" % natoms)
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

    count=0
    linesperframe=natoms+2
    print("Writing out files")
    os.makedirs("Farm",exist_ok=True)
    with open(params["filename"]) as f:
        while(count < params["stop"]/params["framefreq"]):
            os.makedirs("Farm/"+str(count),exist_ok=True)
            start = count*params["framefreq"]
            end = (count+1)*params["framefreq"]
            print("Writing %d directory of %d" % (count,int(params["stop"]/params["framefreq"])))
            g=open("Farm/"+str(count)+"/traj.xyz","w")
            for frame in range(params["framefreq"]):
                for l in range(linesperframe):
                    g.write(f.readline())
            g.close()
            np.savetxt("Farm/"+str(count)+"/e_init.out", np.c_[energy["e"][start:end]])
            np.savetxt("Farm/"+str(count)+"/L.dat", np.c_[Lval[start:end]])
            np.savetxt("Farm/"+str(count)+"/lj_init.out", np.c_[energy["lj"][start:end]])
            np.savetxt("Farm/"+str(count)+"/vol_init.out", np.c_[energy["vol"][start:end]])
            np.savetxt("Farm/"+str(count)+"/ke_init.out", np.c_[energy["ke"][start:end]])
            count+=1
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
    params["h_o"],params["h_h"]=init_rdf(params,len(o)),init_rdf(params,len(h))
    data={ "Roo":[], "Roh":[], "Ctheta":[], "Asphere":[],'Groo':[], 'Groh':[], 'Grhh':[]}
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
            roo,roh,ctheta,asp,goo,goh,ghh=do_analysis(params,frame)
            # The following section increments storage vectors
            data["Roo"].append(roo)
            data["Roh"].append(roh)
            data["Ctheta"].append(ctheta*180.0/np.pi)
            data["Asphere"].append(asp)
            data["Groo"].append(goo)
            data["Groh"].append(goh)
            data["Grhh"].append(ghh)
            # End storage Vector section
            framecount+=1
            end = time.time()
            if (framecount % 10 == 0): print("frame: %s \ntime_per_frame: %s seconds\ntotal_time: %s seconds" % (framecount,(end-start)/framecount, end-start))
            # Writes a restart file
            if (framecount % params["Rfreq"] == 0):
                g = open("restart_distribution.data.pkl","wb")
                pickle.dump(data,g)
                g.close()
                h = open("restart_distribution.energy.pkl","wb")
                pickle.dump(energy,h)
                h.close()
            if (framecount % params["stop"] == 0): break
        manipulate_data(params,data,energy)

    return

    

# Input Arguments
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
parser.add_argument('-restart', default=0, help='Boolean value, 2 to start a farm job, 1 to read in pkl files, 0 to not')
parser.add_argument('-restno', default=10, help='Integer value, number of subdirectories to read restarts')
parser.add_argument('-rest_freq', default=1000, help='Integer value, number of configs per restart file')
parser.add_argument('-frame_freq', default=100, help='Integer value, number of configs per subdir')
parser.add_argument('-dr', default=0.1, help='Float value, bin width for rdf')
parser.add_argument('-rmax', default=10.0, help='Float value, end of rdf')
parser.add_argument('-RDF', default=0, help='Integer value, 1 is on, 0 is off')

args = parser.parse_args()

# Constants
kb=0.0019872041

# Dictionary that holds all the information about the simulation run.
inputparams={"filename":str(args.f), "stop":int(args.nconfigs), "htype":int(args.hatom), "otype":int(args.oatom), "ovtype":str(args.order), "nblocks":int(args.nblocks), "T":float(args.T), "P":float(args.P),"pre":str(args.prepend), "aspon":int(args.asphere),"R":int(args.restart), "numR":int(args.restno),'Rfreq':int(args.rest_freq),'framefreq':int(args.frame_freq),'dr':float(args.dr), 'rmax':float(args.rmax),'rdfon':int(args.RDF)}
inputparams["rdfbins"]=int(inputparams["rmax"]/inputparams["dr"])

print("Welcome to the Distribution Predictor!")
if(inputparams["aspon"]==1): print("Note: Asphericity calculation is on, calcualtion will be much slower.")
if(inputparams["rdfon"]==1): print("Note: RDF calculation is on, calcualtion will be much slower.")


if (inputparams["R"] == 0):
    print("Beginning to read trajectory")
    read_traj(inputparams)
elif (inputparams["R"] == 1):
    print("Restart has been selected")
    call_calculation(inputparams)
elif (inputparams["R"] == 2):
    print("Restart farming is initiated")
    farm_frames(inputparams)
else:
    print("Restart option must be 1 or 0, nothing else will suffice.")
    print("Please try again")
