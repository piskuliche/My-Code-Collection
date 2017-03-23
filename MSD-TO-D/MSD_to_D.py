#!/usr/bin/env python

#This is a code meant for calculating the value of the diffusion coefficient based on the slope of the linear regime of the mean squared displacement
#This code was created by Ezekiel Piskulich at the University of Kansas, 2017. All rights reserved.


import numpy as np
import argparse
from argparse import RawTextHelpFormatter,SUPPRESS

parser = argparse.ArgumentParser(description='''Calculates the diffusion coefficient.''',formatter_class=RawTextHelpFormatter) 

parser.add_argument('-f',help='filename')
args = parser.parse_args()
filename=args.f
filenameout=filename+"dif.dat"
data1,data2,data3=np.genfromtxt(filename,dtype=float,usecols=(0,1,2), autostrip=True,unpack=True)
out=open(filenameout,'w')
for i in range(len(data1)):
    out.write("%f %f %f\n" % (float(data1[i]),float(data2[i]/60*10**5),float(data3[i]/60*10**5)))
    
