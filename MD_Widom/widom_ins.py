#!/usr/bin/env python
"""
This is a python program that can be called to generate trial configurations for lammps
"""
import numpy as np
import sys
import conv_connectivity as cc


def read_xyz(connectfile):
    # Reads in coordinates from a file
    xyz_fname = connectfile.replace(".connect", ".xyz")
    x, y, z = np.genfromtxt(xyz_fname, usecols=(0,1,2), skip_header = 2, unpack=True)
    return x, y, z

def calculate_com(r, M):
    com_r = [0,0,0]
    natms = len(x)
    for i in range(3):
        for atm in range(natms):
            com_r[i] += M[i]*r[i][atm]
    com_r[:] = com_r[:]/np.sum(M)
    return com_r

def new_orientation(r, com_r, theta, phi):
    # This function uses rotation matrices to create a new orientation
    

        
        


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Error: improper input arguments")
        print("Usage: python widom_ins.py connectfile")
        sys.exit()
    else:
        connectfile = str(sys.argv[1])
        print("Read in connectivity file")
    natms, M, header, coords, footer = cc.read_connectivity(connectfile)
    print("There are %d atoms" % natms)
    x, y, z = cc.read_xyz(connectfile)
    r = [x,y,z]
    assert len(x) == natms, "Error: Number of atoms different than number of coordinates"
    cc.write_connectivity(x, y, z, header, coords, footer, connectfile)
    print("Molec file has been generated")
