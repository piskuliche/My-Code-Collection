import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

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
            if dist < rmax+dr/2:
                counts[_find_nearest(xrdf,dist)] += 2
                fcount += 1
        h_real = h_real + counts/natoms
    # Calculate ideal value
    h_id = math.pi*natoms/(L[0]*L[1])*(rhi**2-rlo**2)
    rdf = h_real/h_id
    # Finalize rdf normalization
    return xrdf, rdf,h_real
                



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calc RDF')
    parser.add_argument('-fcount',default=10000,type=int,help='Number of frames [default 10000]')
    parser.add_argument('-dim',default=2,type=int,help='Number of dimensions [default 2]')
    parser.add_argument('-dr', default=0.1, type=float, help="Bin thickness [default 0.1]")
    parser.add_argument('-infile', default="DPPC.fluid.lammpsdump", type=str, help="File trajectory name [default DPPC.fluid.lammpsdump]")
    parser.add_argument('-outfile',default="tot_rdf.out",type=str,help='Output file name [default tot_rdf.out]')
    args = parser.parse_args()
    
    infile=args.infile
    outfile=args.outfile
    dr=args.dr
    dim=args.dim
    fcount = args.fcount


    u = mda.Universe("system.data",infile)
    count = 0
    tot_rdf = []
    tot_h = []
    for ts in u.trajectory:
        if count%1000 == 0:
            print("Current count is %d" % count)
        leafs = LeafletFinder(u, 'type 4')
        leaf1,leaf2 = leafs.groups(0), leafs.groups(1) 
        r,rdf,h_real = calc_RDF(u,leaf1,dr=dr)
        tot_rdf.append(rdf)
        tot_h.append(h_real)
        r,rdf,h_real = calc_RDF(u,leaf2,dr=dr)
        tot_rdf.append(rdf)
        tot_h.append(h_real)
        count += 1
        if count == fcount:
            np.savetxt(outfile,np.c_[r,np.average(tot_rdf,axis=0),np.average(tot_h,axis=0)])
            exit()
