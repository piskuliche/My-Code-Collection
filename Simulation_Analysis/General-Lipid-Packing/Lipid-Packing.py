#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
import argparse

from MDAnalysis.analysis.leaflet import LeafletFinder
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch as Neighbors
from numpy.linalg import norm
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d



def _pbc_dist(box,r1,r2):
    """
    This is a general distance calculator that calculates minimum
    image distances. 

    Input: 
    -box (6 element array in the MDAnalysis style)
    -r1  (3 element coordinates)
    -r2  (3 element coordinates)

    Output:
    -dr  (3 element distance unit vector)
    """
    L = np.array(box[:3])
    dr = r1-r2
    dr = dr + L*np.array(list(map(int,dr/L)))
    dr /= norm(dr,2)
    return dr

def Select_Resids(u,Iargs):
    """
    Function to select residues of interest

    Note - this may need to be modified for particular purposes depending on atom selections

    Input:
    -u: MDAnalysis universe
    -Iargs: input args

    Output: MDAnalysis selections
    """
    building_blocks = {}
    with open(Iargs.keys,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            if len(line) > 1:
                building_blocks[line[0]]=list(map(int, line[1:]))
    start = 1
    stop  = building_blocks[Iargs.bb[0]][0]
    if len(Iargs.bb) == 2:
        start = np.sum(building_blocks[Iargs.bb[0]])+1
        stop  = start + building_blocks[Iargs.bb[1]][0]-1
    print("Start: %d Stop: %d" % (start,stop))
    PC_sel = u.select_atoms("resid %d to %d" % (start,stop))
    PO4_sel = u.select_atoms("resid %d to %d and type 6" % (start,stop))
    return PC_sel, PO4_sel

def Get_Leaflets(u,PO4_sel):
    # Find leaflets
    L = LeafletFinder(u,PO4_sel,pbc=True,sparse=False)
    leaf_sizes = L.sizes()
    # Store Sizes
    n_per_leaf, key_index = [],[]
    for key in leaf_sizes:
        n_per_leaf.append(leaf_sizes[key])
        key_index.append(key)
    # Sort
    n_leaf_srt, n_key_srt = zip(*sorted(zip(n_per_leaf,key_index)))
    # Group Selections
    in_PO4_sel  = L.groups(n_key_srt[-2:][0])
    out_PO4_sel = L.groups(n_key_srt[-2:][1])
    in_sel      = u.select_atoms("same resid as group leaf", leaf=in_PO4_sel)
    out_sel     = u.select_atoms("same resid as group leaf", leaf=out_PO4_sel)
    
    return in_sel, out_sel, n_leaf_srt

def voronoi_volumes(v):
    """
    This is a function pulled from here:
    https://stackoverflow.com/a/54265278/11361102

    It sorts the voronoi volumes by the point_region
    and returns the volume
    """
    from scipy.spatial import ConvexHull
    vol = np.zeros(v.npoints)
    for i, reg_num in enumerate(v.point_region):
        indices = v.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(v.vertices[indices]).volume
    return vol

def Lipid_Pointer_and_Neighbors(u,selection,Iargs):
    """
    This is the meat of the calculation. It loops over lipids and calculates the 
    area for each lipid using voronoi tesselation.
    
    """
    def Project_Points(points3d, origin, normal, u, v):
        """
        Simple function for projecting points from 3d orthogonally onto a 2d plane
        It reports them in planar coordinates
        """
        points2d=[]
        for point in points3d:
            d = np.dot((point-origin),normal)
            proj3d = point-d*normal
            projx = np.dot(u,proj3d-origin)
            projy = np.dot(v,proj3d-origin)
            projz = np.dot(normal,proj3d-origin)
            points2d.append([projx,projy])
        return np.array(points2d)
    box = u.dimensions
    sel = selection.split('residue')
    vectors=[]
    #fig = plt.figure()
    
    #ax = plt.axes(projection='3d')
    vols = []
    for lipid in sel:

        pointer=_pbc_dist(box,lipid.positions[1],lipid.positions[-1])
        vectors.append(pointer)
        neigh = Neighbors(selection)

        patom = lipid.select_atoms("%s" % Iargs.pselect)
        sval = neigh.search(patom,15.0,level='R')
        sel_P =  sval.atoms.select_atoms("%s" % Iargs.pselect)

        svd = np.linalg.svd(sel_P.positions.T - np.mean(sel_P.positions.T, axis=1, keepdims=True))
        left=svd[0]
        normal_vec, pvec1, pvec2 = left[:,2], left[:,0], left[:,1]
        pos = patom.positions[0]

        #ax.scatter3D(sel_P.positions.T[0],sel_P.positions.T[1],sel_P.positions.T[2])
        #ax.quiver(pos[0],pos[1],pos[2],normal_vec[0],normal_vec[1],normal_vec[2],length=5) 
        # Do the Projection
        proj2d = Project_Points(sel_P.positions,pos,normal_vec,pvec1,pvec2)

        vor = Voronoi(proj2d)
        vol = voronoi_volumes(vor)
        #voronoi_plot_2d(vor)
        #plt.show()

        voltmp = vol[np.where(proj2d.T[0]==0)][0]
        if voltmp != np.inf:
            vols.append(voltmp)
    #print(np.average(vols))  
    #plt.show()
    
    return vols



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-data',    default="equil.data",       type=str, help='Data file name [default=equil.data]')
    parser.add_argument('-start',   default=0,                  type=int, help='Starting file')
    parser.add_argument('-stop',    default=11,                 type=int, help='Stopping file')
    parser.add_argument('-pre',     default='coords.lmpstrj-',  type=str, help='File prefix')
    parser.add_argument('-keys',    default='atom.keys',        type=str, help='File with atom numbers inside')
    parser.add_argument('-bb',      nargs='+',                  type=str, help='Names of building blocks')
    parser.add_argument('-skip',    default=20,                 type=int, help='Number of frames to skip between frames')
    parser.add_argument('-histbins',default=50,                 type=int, help='Number of bins for histogramming')
    parser.add_argument('-subdir',  default="y",                type=str, help = 'Dumpfiles in numbered subdirs y/n')
    parser.add_argument('-pselect', default='type 6',           type=str, help = 'Selection language for P atom')
    Iargs = parser.parse_args()

    # Load Trajectory Files
    trjfiles = []
    for i in range(Iargs.start,Iargs.stop+1):
        tmp_tuple = None
        if Iargs.subdir == "y":
            tmp_tuple = (str(i) + "/"+Iargs.pre + str(i),'LAMMPSDUMP')
        else if Iargs.subdir == "n":
            tmp_tuple = (Iargs.pre + str(i),'LAMMPSDUMP')
        else:
            exit("Error: improper choice of subdir")
        trjfiles.append(tmp_tuple)

    # Load MDAnalysis Universe
    u = mda.Universe(Iargs.data,trjfiles,dt=200.0)

    print("Data has been read")
    
    # Choose residues
    PC_sel,PO4_sel = Select_Resids(u,Iargs)
    vol_in, vol_out = [],[]
    N_in, N_out = [],[]
    # Loop over Frames
    for ts in u.trajectory[::Iargs.skip]:
        print("Frame: %d" % ts.frame)
        # Get Leaflets
        in_sel, out_sel, nsrt = Get_Leaflets(u,PO4_sel)

        N_in.append(nsrt[0])
        N_out.append(nsrt[1])

        # Find voluems
        tmp_vols_in=Lipid_Pointer_and_Neighbors(u,in_sel)
        tmp_vols_out=Lipid_Pointer_and_Neighbors(u,out_sel)
        vol_in  = np.append(vol_in , tmp_vols_in)
        vol_out = np.append(vol_out, tmp_vols_out)

    # Do plotting but also save files with volumes and leaflets.
    plt.hist(vol_in,bins=Iargs.histbins,range=(30,200), density=True)
    plt.show()
    print(np.average(vol_in)/100)
    np.savetxt("inner_lipid_volumes.dat",np.c_[vol_in])
    plt.hist(vol_out,bins=Iargs.histbins,range=(30,200), density=True)
    plt.show()
    print(np.average(vol_out)/100)
    np.savetxt("outer_lipid_volumes.dat",np.c_[vol_out])
    np.savetxt("N_per_leaf.dat",np.c_[N_in, N_out])
    
    









