#!/usr/bin/env python
"""
This is a python program that can be called to generate trial configurations for lammps
"""
import numpy as np
import sys
import conv_connectivity as cc
import random

def init_random():
    random.seed(a=None)
    return

def read_xyz(connectfile):
    # Reads in coordinates from a file
    xyz_fname = connectfile.replace(".connect", ".xyz")
    x, y, z = np.genfromtxt(xyz_fname, usecols=(0,1,2), skip_header = 2, unpack=True)
    return x, y, z

def calculate_com(r, M):
    com_r = [0,0,0]
    natms = len(r)
    for i in range(3):
        for atm in range(natms):
            com_r[i] += M[i]*r[atm][i]
    com_r[:] = com_r[:]/np.sum(M)
    return com_r



def rand_quaternion():
    """This implements the algorithm of Marsaglia 1972 and Vesely 1982 to calculate a random unit quaternion"""
    # initialize random numbers
    rval =[0.0,0.0,0.0,0.0]
    s1 = 5
    # Pick first two random numbers such that s1 = r1^2 + r2^2 < 1
    while s1 > 1:
        rval[0]=2*random.random()-1
        rval[1]=2*random.random()-1
        s1 = rval[0]**2. + rval[1]**2.
    s2 = 5
    # Pick second two random numbers such that s1 = r1^2 + r2^2 < 1
    while s2 > 1:
        rval[2]=2*random.random()-1
        rval[3]=2*random.random()-1
        s2 = rval[2]**2. + rval[3]**2.
    # Define quaternion
    s_w = np.sqrt((1-s1)/s2)
    a = [rval[0],rval[1],rval[2]*s_w, rval[3]*s_w]
    asq = a[0]**2.+a[1]**2.+a[2]**2.+a[3]**2.
    assert  round(asq,5) == 1, ("A is not a unit vector %s" % asq)
    return a
    
def quaternion_to_rotmatrix(a):
    """ This generates the rotation matrix from the quaternion of interest"""
    A11 = a[0]**2.+a[1]**2.-a[2]**2.-a[3]**2.
    A22 = a[0]**2.-a[1]**2.+a[2]**2.-a[3]**2.
    A33 = a[0]**2.-a[1]**2.-a[2]**2.+a[3]**2.
    A12 = 2*(a[1]*a[2]+a[0]*a[3])
    A13 = 2*(a[1]*a[3]-a[0]*a[2])
    A21 = 2*(a[1]*a[2]-a[0]*a[3])
    A23 = 2*(a[2]*a[3]+a[0]*a[1])
    A31 = 2*(a[1]*a[3]+a[0]*a[2])
    A32 = 2*(a[2]*a[3]-a[0]*a[1])

    A = [[A11,A12,A13],[A21,A22,A23],[A31,A32,A33]]
    return A

def rand_rot(r, com_r):
    """Chooses a random quaternion and then applies it to the molecule"""
    a = rand_quaternion()
    A = quaternion_to_rotmatrix(a)
    newr = []
    for r_atom in r:
        dr = []
        for i in range(len(r_atom)):
            dr.append(r_atom[i]-com_r[i])
        tmpr = np.matmul(A, dr) + com_r
        newr.append(tmpr)
    return newr

def print_rot(r,newr):
    """prints the old and the new rotations"""
    f=open("orig.xyz",'w')
    f.write("%d\n" % (len(r)+1))
    f.write("original\n")
    f.write("O %.5f %.5f %.5f\n" %(0.0, 0.0, 0.0))
    for r_atom in r:
        f.write("H %.5f %.5f %.5f\n" %(r_atom[0], r_atom[1], r_atom[2]))
    g=open("new.xyz",'w')
    g.write("%d\n" % (len(newr)+1))
    g.write("new\n")
    g.write("O %.5f %.5f %.5f\n" %(0.0, 0.0, 0.0))
    for newr_atom in newr:
        g.write("H %.5f %.5f %.5f\n" %(newr_atom[0], newr_atom[1], newr_atom[2]))
    f.close()
    g.close()
    return

def rand_trans(r,com_r,L):
    rval = [0,0,0]
    newcom_r = []
    r_trans = []

    rval[0] = random.random()
    rval[1] = random.random()
    rval[2] = random.random()

    for i in range(len(com_r)):
        newcom_r.append(rval[i]*L)
    
    for atom in range(len(r)):
        r_shift = []
        for i in range(3):
            tmpr = r[atom][i] - com_r[i] + newcom_r[i]
            r_shift.append(tmpr)
        r_trans.append(r_shift)
    return r_trans, newcom_r



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
    r=np.transpose([x,y,z])
    com_r =calculate_com(r, M)
    print(com_r)
    assert len(x) == natms, "Error: Number of atoms different than number of coordinates"
    cc.write_connectivity(x, y, z, header, coords, footer, connectfile)
    print("Molec file has been generated")
    init_random()
    r_rotated = rand_rot(r,com_r)
    r_translated, com_r = rand_trans(r_rotated,com_r,L)

