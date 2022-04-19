#!/usr/bin/env python

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
from scipy.special import sph_harm

def calc_Q6(ts, leaf,Nb=5):
    """
    This function calculates the Q6 order parameter for lipids in a leaflet. 

    <Q_6> = <[4pi/13 sum_{-6}^{6}(1/Nb sum_j^Nb Y_6m(theta,phi))^2]^(1/2)>

    Note theta = pi/2
    phi = angle between arbitrary coordinate system and neighbor. 
    Nb = number nearest neighbors
    
    """
    def _pbc_dist(r1,r2,L,dim=2):
        dr = np.zeros(dim)
        drsq = 0
        for i in range(dim):
            dr[i] = r1[i]-r2[i]
            dr[i] -= L[i]*round(dr[i]/L[i])
            drsq += dr[i]**2
        dist = np.sqrt(drsq)
        return dist, dr

    L = ts.dimensions[:2]
    const = 4*math.pi/13
    mval = 6
    natoms=np.shape(leaf.positions)[0]

    # Find the nearest neighbors
    distmat = np.ones((natoms,Nb))*5e5 
    ref_vec = np.array([1,0])
    angmat = np.zeros((natoms,Nb))
    for i in range(natoms):
        ri = leaf.positions[i]
        for j in range(natoms):
            rj = leaf.positions[j]
            if i != j:
                rij_dist, dr  = _pbc_dist(ri,rj,L)
                if rij_dist < np.max(distmat[i]):
                    idx = np.argmax(distmat[i])
                    distmat[i][idx] = rij_dist
                    urj = np.divide(dr,rij_dist)
                    ang = np.arccos(np.dot(ref_vec,urj))
                    angmat[i][idx]=ang
    Q6 = 0
    for i in range(natoms):
        m_sum_values = []
        for m in range(-mval,mval+1):
            jsum = 0
            for j in range(0,Nb):
                jsum += sph_harm(m,6,angmat[i,j],math.pi/2)
            m_sum_values.append(np.abs((jsum/Nb))**2)
        Q6 += np.sqrt(np.sum(m_sum_values)*const)
    Q6 = Q6/natoms

    return Q6
                



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calc RDF')
    parser.add_argument('-Nb', default=5, type=int, help="Number of nearest neighbors")
    parser.add_argument('-infile', default="DPPC.fluid.lammpsdump", type=str, help="File trajectory name")
    parser.add_argument('-outfile',default="Q6.out",type=str,help='Output file name')
    parser.add_argument('-fcount',default=10000,type=int,help='Number of configurations')
    args = parser.parse_args()
    
    infile=args.infile
    outfile=args.outfile
    Nb=args.Nb
    fcount=args.fcount

    u = mda.Universe("system.data",infile)
    count = 0
    Qsix=[]
    for ts in u.trajectory:
        if count%1000 == 0:
            print("Current count is %d" % count)
        leafs = LeafletFinder(u, 'type 4')
        leaf1,leaf2 = leafs.groups(0), leafs.groups(1) 
        Qsix.append(calc_Q6(u,leaf1,Nb=Nb))
        Qsix.append(calc_Q6(u,leaf2,Nb=Nb))
        count += 1
        f=open(outfile,'w')
        if count == fcount:
            f.write("%s %10.5f\n" % (infile,np.average(Qsix)))
            f.close()
            exit()
