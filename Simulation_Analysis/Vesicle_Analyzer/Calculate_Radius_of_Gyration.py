#!/usr/bin/env python
import MDAnalysis as mda
import numpy as np
import argparse
from MDAnalysis.analysis.base import AnalysisBase,AnalysisFromFunction,analysis_class
from MDAnalysis.analysis.leaflet import LeafletFinder
from scipy.spatial import SphericalVoronoi

parser = argparse.ArgumentParser()
parser.add_argument('-data', default="equil.data", type=str, help='Data file name [default=equil.data]')
parser.add_argument('-start', default=0, type=int, help='Starting file')
parser.add_argument('-stop', default=11, type=int, help='Stopping file')
parser.add_argument('-pre', default='coords.lmpstrj-', type=str, help='File prefix')
parser.add_argument('-keys',default='atom.keys', type=str, help='File with atom numbers inside')
parser.add_argument('-bb',nargs='+', type=str, help='Names of building blocks')
args= parser.parse_args()

datafile = args.data
trjfile = args.pre
num_start = args.start
num_stop  = args.stop
keyfile = args.keys
bb_choice = args.bb

def radgyr(atomgroup, masses, total_mass=None):
    """
    This is a function to calculate Rg, by default it is not actually used,
    mostly because PBC corrections need to be included first. 
    """
    # coordinates change for each frame
    coordinates = atomgroup.positions
    center_of_mass = atomgroup.center_of_mass()

    # get squared distance from center
    ri_sq = (coordinates-center_of_mass)**2
    # sum the unweighted positions
    sq = np.sum(ri_sq, axis=1)
    sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
    sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
    sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y

    # make into array
    sq_rs = np.array([sq, sq_x, sq_y, sq_z])

    # weight positions
    rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
    # square root and return
    return np.sqrt(rog_sq)

# Read building block info from atom.keys file
building_blocks={}
with open(keyfile,'r') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split()
        if len(line) > 1:
            building_blocks[line[0]]=list(map(int, line[1:]))


# Make list of trajectory file names
trjfiles=[]
for i in range(num_start, num_stop+1):
    tmp_tuple = (trjfile+str(i),'LAMMPSDUMP')
    trjfiles.append(tmp_tuple)

u = mda.Universe(datafile,trjfiles,dt=200.0)

print("Data has been read")
select_PS = u.select_atoms("resid 1")
# Need to figure out initial resid
print(bb_choice)
start_residue = np.sum(building_blocks[bb_choice[0]])+1
stop_residue = start_residue + building_blocks[bb_choice[1]][0]-1
# Final Selections
select_PC = u.select_atoms("resid %d to %d" % (start_residue,stop_residue))
select_all = u.select_atoms("resid 1 %d-%d" % (start_residue,stop_residue))
select_PC_PO4 = u.select_atoms("resid %d to %d and type 8" % (start_residue,stop_residue))

print("PS has %d atoms" % select_PS.n_atoms)
print("PC has %d atoms" % select_PC.n_atoms)
print("PS+PC has %d atoms" % select_all.n_atoms)

timesteps, PC_Rgs, PS_Rgs, All_Rgs = [],[],[],[]
lfone,lftwo,num_leafs,lfother = [],[],[],[]
print("Beginning loop")
for ts in u.trajectory[::100]:
    print(ts.time)
    # Find Leaflets
    L = LeafletFinder(u,select_PC_PO4,pbc=True,sparse=False)
    # Count total num leaflets
    leaf_sizes = L.sizes()
    print(leaf_sizes)
    # Create list of leaflets,keys
    n_per_leaf,key_index = [],[]
    for key in leaf_sizes:
        n_per_leaf.append(leaf_sizes[key])
        key_index.append(key)
    # Sort leaflets, keys
    n_leaf_srt, n_key_srt = zip(*sorted(zip(n_per_leaf,key_index)))
    lf1 = L.groups(key_index[-2:][0])
    lf2 = L.groups(key_index[-2:][1])
    L.write_selection("test.vmd")
    tmp_n_lf1, tmp_n_lf2 = n_leaf_srt[-2:][0],n_leaf_srt[-2:][1]
    tmp_other = np.sum(n_per_leaf) - tmp_n_lf1 - tmp_n_lf2
    lfone.append(tmp_n_lf1)
    lftwo.append(tmp_n_lf2)
    lfother.append(tmp_other)
    print(lf1)
    print(np.shape(lf1.positions))
    com=L.groups(0).center_of_mass()
    rad=select_PC_PO4.radius_of_gyration(pbc=True)
    #SphericalVoronoi(L.groups(0).positions,radius=rad,center=com)
    num_leafs.append(len(n_per_leaf))
    timesteps.append(ts.time)
    PC_Rgs.append(select_PC.radius_of_gyration(pbc=False))
    PS_Rgs.append(select_PS.radius_of_gyration(pbc=False))
    All_Rgs.append(select_all.radius_of_gyration(pbc=False))

np.savetxt("All-Rgs.dat",np.c_[timesteps,PC_Rgs,PS_Rgs,All_Rgs])
np.savetxt("All-lfs.dat", np.c_[timesteps,lfone,lftwo,lfother,num_leafs])
