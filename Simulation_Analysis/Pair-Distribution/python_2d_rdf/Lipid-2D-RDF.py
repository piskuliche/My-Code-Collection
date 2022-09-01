#!/usr/bin/env python
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import matplotlib.pyplot as plt
import math,os
import argparse
import pickle

def calc_RDF(ts,leaf,dr=0.1,rmax=18,dim=2):
    """
    This is a function that takes in an MDAnalysis timestep, and leaflet AtomGroup object, and returns
    the radial distribution function calculated from the frame. 

    Inputs:
        -timestep object
        -leaflet AtomGroup object
        -dr is the shell thickness for the rdf calculation
        -rmax is the maximum distance to run the rdf calculation over. 

    Outputs:
        -xrdf: array of shell midpoints
        -rdf: array of rdf values
        -h_real: array of individual bin counts.
    """
    def _pbc_dist(r1,r2,L,dim=2):
        """
        Simple function that takes two vectors and calculates
        the minimum image distance between the two.

        Inputs:
            -r1, r2: positions arrays for each of the vectors shape(dim)
            -L: box length array of shape(dim)
            -dim: # of dimensions
        Output:
            -dist: minimum image distance.
        """
        dr = np.zeros(dim)
        drsq = 0
        # loop over dimensions
        for i in range(dim):
            dr[i] = r1[i]-r2[i]
            dr[i] -= L[i]*round(dr[i]/L[i])
            drsq += dr[i]**2
        # output sqrt distance
        dist = np.sqrt(drsq)
        return dist
    def _find_nearest(array, value):
        """
        Finds the nearest location in array to value.
        returns: index of location.
        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    # Calculate # bins
    nbins = int(rmax/dr)
    # Select # atoms
    natoms=np.shape(leaf.positions)[0]
    # Box dimensions
    L=ts.dimensions[:dim]

    # Set up Arrays
    h_real = np.zeros(nbins)
    xrdf = np.arange(0,rmax,dr)+(dr/2)
    rlo = np.arange(0,rmax,dr)
    rhi = rlo+dr
    nfinal = natoms
    
    # Loop over atoms
    for i in range(natoms):
        counts = np.zeros(nbins)
        fcount = 0
        # Loop over i+1, natoms
        for j in range(i+1,natoms):
            dist=_pbc_dist(leaf.positions[i],leaf.positions[j],L)
            if dist < rmax:
                counts[_find_nearest(xrdf,dist)] += 2
                fcount += 1
        h_real = h_real + counts/natoms
    # Calculate ideal value
    h_id = math.pi*natoms/(L[0]*L[1])*(rhi**2-rlo**2)
    rdf = h_real/h_id
    # Finalize rdf normalization
    return xrdf, rdf,h_real
               
def Weight_RDFs(r,rdfs,efile):
    fnames, etypes = np.genfromtxt(efile,usecols=(0,1),dtype=str,unpack=True)
    energies={}
    nrdf=np.shape(rdfs)[0]
    for i in range(len(etypes)):
        energies[etypes[i]] = np.genfromtxt(fnames[i],usecols=0)
        nskip = int(len(energies[etypes[i]])/nrdf)
        energies[etypes[i]]=energies[etypes[i]][::nskip]
        de = (energies[etypes[i]]-np.average(energies[etypes[i]]))[:-1]
        erdf = -np.multiply(de[:,None],rdfs)
        np.savetxt("%s_%s_rdf.out"%(label,etypes[i]),np.c_[r,np.average(erdf,axis=0)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calc RDF')
    parser.add_argument('-fcount',default=10000,type=int,help='Number of frames per segment [default 1000]')
    parser.add_argument('-dim',default=2,type=int,help='Number of dimensions [default 2]')
    parser.add_argument('-dr', default=0.1, type=float, help="Bin thickness [default 0.1]")
    parser.add_argument('-data', default="system.data", type=str, help="Data file [default system.data]")
    parser.add_argument('-infile', default="DPPC.fluid.lammpsdump", type=str, help="File trajectory name [default DPPC.fluid.lammpsdump]")
    parser.add_argument('-outfile',default="tot_rdf.out",type=str,help='Output file name [default tot_rdf.out]')
    parser.add_argument('-leafselec',default='type 4',type=str, help='Leaflet selection text [default "type 4"]')
    parser.add_argument('-segid',default=0,type=int,help='Which segment to do?')
    parser.add_argument('-rmax',default=18,type=float,help='Maximum distance in r')
    parser.add_argument('-op',default=0,type=int,help='Which option? [0] calc rdf [1] sum rdf')
    parser.add_argument('-ew',default=0,type=int,help='[0] No EWeighting, [1] Eweighting')
    parser.add_argument('-efile',default="flucts.inp", help='Name of file with energy types')
    parser.add_argument('-label',default="PO4",help="Label to identify files")
    args = parser.parse_args()
    
    infile=args.infile
    outfile=args.outfile
    dr=args.dr
    dim=args.dim
    fcount = args.fcount
    leafselec = args.leafselec
    datafile = args.data
    segid = args.segid
    rmax = args.rmax
    op = args.op
    label = args.label
    ew = args.ew
    efile = args.efile
    path="%s_rdfs"%label
    if op == 0:
        ew = 0
        efile = None
        os.makedirs("%s_rdfs"%label,exist_ok=True)

    if op == 0:
        u = mda.Universe(datafile,infile)
        count = 0
        tot_rdf = []
        tot_h = []
        segstart,segstop = segid*fcount, (segid+1)*fcount
        # Loop over trajectory
        for ts in u.trajectory[segstart:segstop]:
            print(ts.frame)
            leafs = LeafletFinder(u, leafselec,pbc=True)
            leaf1,leaf2 = leafs.groups(0), leafs.groups(1)
            print("leaf1")
            r,rdf1,h_real1 = calc_RDF(u,leaf1,dr=dr,rmax=rmax)
            print("leaf2")
            r,rdf2,h_real2 = calc_RDF(u,leaf2,dr=dr,rmax=rmax)
            tot_rdf.append((rdf1+rdf2)/2)
            tot_h.append((h_real1+h_real2)/2)
            pickle.dump(tot_rdf,open(path+"/%s_rdf_%d.pckl"%(label,segid),'wb'))
            pickle.dump(tot_h,open(path+"/%s_h_%d.pckl"%(label,segid),'wb'))
            if segid == 0:
                pickle.dump(r,open(path+"/%s_r.pckl"%label,'wb'))
    else:
        tot_rdf = []
        tot_h = []
        r = pickle.load(open(path+"/%s_r.pckl"%label,'rb'))
        for seg in range(segid):
            rdfs=pickle.load(open(path+"/%s_rdf_%d.pckl"%(label,seg),'rb'))
            hs=pickle.load(open(path+"/%s_h_%d.pckl"%(label,seg),'rb'))
            print(seg,np.shape(rdfs),np.shape(hs))
            for rdf in rdfs:
                tot_rdf.append(rdf)
            for h in hs:
                tot_h.append(h)
        np.savetxt(label+"_"+outfile,np.c_[r,np.average(tot_rdf,axis=0),np.average(tot_h,axis=0)])
        if ew == 1: 
            Weight_RDFs(r,tot_rdf,efile)

