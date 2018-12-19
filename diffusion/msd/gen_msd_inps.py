#!/usr/bin/env python
'''
This is a general python code to create the input for the msd calculation c++ code.
Usage: gen_msd_inps.py datafile
'''
import numpy as np
import sys


# System Arguments
if len(sys.argv) != 2:
    print("Usage: gen_msd_inps.py datafile ")
    exit(0)
datafile = str(sys.argv[1])

option = input("[1] NVT Run [2] NPT Run -> NVT Run\n")
# Read Data File for Dimensions
if option == 1:
    lines=[]
    with open(datafile, 'r') as f:
        for line in f:
            if "xlo" in line:
                lines.append(line.split()[:2])
            elif "ylo" in line:
                lines.append(line.split()[:2])
            elif "zlo" in line:
                lines.append(line.split()[:2])
                break
if option === 2:
    print("Write out the box.info file yourself! Don't be lazy!")
    print("Okay... yeah I am being lazy right now as I write this")

# Generates box.info file with dimensions
boxfile = open("box.info", 'w')
for line in lines:
    boxfile.write("%s %s\n" % (line[0],line[1]))
boxfile.close()

# Calls for user response to create molecule input file for msd code
print("Now generating a molecule input file")
name            = str(raw_input("What is the name of the molecule?\n"))
atoms_per_mol   = input("How many atoms per mol?\n")
msd_index       = input("What msdindex would you like to do (between 1 and %s)\n" % atoms_per_mol)
c2index1        = input("What c2index1 would you like to do (between 1 and %s)\n" % atoms_per_mol)
c2index2        = input("What c2index2 would you like to do (between 1 and %s)\n" % atoms_per_mol)
m = []
for i in range(atoms_per_mol):
    m.append(float(raw_input("What is the mass of atom %s\n" % (i+1))))

# Writes out to a file.
molfile = open(name+".txt", 'w')
molfile.write("%s\n" % atoms_per_mol)
molfile.write("%s\n" % msd_index)
molfile.write("%s %s\n" % (c2index1,c2index2))
for i in range(atoms_per_mol):
    molfile.write("%s\n" % m[i])
