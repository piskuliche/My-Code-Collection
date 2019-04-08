#!/usr/bin/env python

import numpy as np
import argparse
from argparse import RawTextHelpFormatter,SUPPRESS

parser = argparse.ArgumentParser(description='''Finds minimum value of a data file''',formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', help ="File Name")

args = parser.parse_args()
filename=args.f
file1 = open(args.f)
file1.readline()
line1=file1.readline()
line1 = line1.split()
file1.close()
print line1, '\n'
print 'length: ',len(line1), '\n'        
data = np.loadtxt(filename, skiprows=3)
for i in range(len(line1)-1):
    if line1[i+1] == 'Volume':
        y_data = data[:,i]
        print min((y_data ** (1./3))/2.) 
