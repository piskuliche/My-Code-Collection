#!/usr/bin/env python
"""
This is a python module for converting lammps connectivity files
It reads the coordinateless sample connectivity file that is written by my general build code 
    and adds in the coordinates. 


"""

import numpy as np
import sys


def read_connectivity(connectfile):
    # Reads in the barebones connectivity file
    header_lines = []
    coords_lines = []
    footer_lines = []
    mass         = []
    typ          = []
    startflag, endflag = 0, 0
    mstart, mend, mcount = 0, 0, 0
    tstart, tend, tcount = 0, 0 ,0
    natms = 0
    with open(connectfile) as cnf:
        for line in cnf:
            if "Coords" in line:
                startflag = 1
            elif "Types" in line:
                endflag   = 1
                tstart    = 1
            elif "Masses" in line:
                mstart = 1
            elif "Bonds" in line:
                mend   = 1
            elif "Charges" in line:
                tend = 1
            if startflag == 0:
                header_lines.append(line)
                if "atoms" in line:
                    natms = int(line.split()[0])
            if startflag == 1 and endflag == 0:
                coords_lines.append(line)
            if endflag == 1:
                assert startflag != 0, "Error: Coords flag missing"
                footer_lines.append(line)
            if mstart == 1 and mend == 0:
                if mcount >= 2 and mcount < 2 + natms:
                    mass.append(float(line.rstrip().split()[1]))
                mcount += 1
            if tstart == 1 and tend == 0:
                if tcount >= 2 and tcount < 2 + natms:
                    typ.append(int(line.rstrip().split()[1]))
                tcount += 1
            
    return natms, mass,typ, header_lines, coords_lines, footer_lines

def read_xyz(connectfile):
    # Reads in coordinates from a file
    xyz_fname = connectfile.replace(".connect", ".newxyz")
    x, y, z = np.genfromtxt(xyz_fname, usecols=(0,1,2), unpack=True)
    return x, y, z

def write_connectivity(x, y, z, header_lines, coords_lines, footer_lines, connectfile):
    # Writes out the new connectivity as a molec file
    assert len(x) == len(coords_lines)-3, "Error: Number of coordinates (%d) different than known connectivity lines (%d)" % (len(x), len(coords_lines)-3)
    molec_fname = connectfile.replace(".connect", ".molec")
    mnf = open(molec_fname, 'w')
    for line in header_lines:
        mnf.write(line)
    for l in range(len(coords_lines)):
        # Note with this condition we don't need to worry about the last blank line
        # due to the fact that it is the same as coord_lines[1]
        line = coords_lines[l]
        if line in coords_lines[:2]:
            mnf.write(line)
        else:
            mnf.write("%s %.5f %.5f %.5f\n" % (line.rstrip(), x[l-2], y[l-2], z[l-2]))
    for line in footer_lines:
        mnf.write(line)
    mnf.close()

        

    

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Error: improper input arguments")
        print("Usage: python ins.py connectfile")
        sys.exit()
    else:
        connectfile = str(sys.argv[1])
        print("Read in connectivity file")
    natms, M,  header, coords, footer = read_connectivity(connectfile)
    print("There are %d atoms" % natms)
    x, y, z = read_xyz(connectfile)
    assert len(x) == natms, "Error: Number of atoms different than number of coordinates"
    write_connectivity(x, y, z, header, coords, footer, connectfile)
    print("Molec file has been generated")
    

   
