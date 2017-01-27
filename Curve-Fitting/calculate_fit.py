#!/usr/bin/python


#This code is a generalized fitting code which should be able to fit most space or tab delimited data. 
#
#Usage: python curve_fit.py -f filename -c1 0 -c1 1 -s1 0 -nl 100 -type exp3
#Copyright Zeke Piskulich, University of Kansas, 2017


import math, sys
import numpy as np
import argparse
from argparse import RawTestHelpFormatter,SUPPRESS
from scipy.optimize import curve_fit


#Argument parser section
parser = argparse.ArgumentParser(description='''Calculates the fit of a general data set''', formatter_class=RawTextHelpFormatter)


parser.add_argument('-f',    help="Data file name")
parser.add_argument('-c1',   help="Numeric Index of Column 0,1,2..N")
parser.add_argument('-c2',   help="Numeric Index of Column 0,1,2..N")
parser.add_argument('-s',    help="Number of lines to skip at top of file")
parser.add_argument('-nl',   help="Number of lines to fit")
parser.add_argument('-type', help="Type of fit to calculate")
parser.add_argument('-i',    help="Extra identifier to append to file")

#Parses the arguments and assigns them to variables
args = parser.parse_args()
file_d=args.f
col1=int(args.c1)
col2=int(args.c2)
skiptop=int(args.s)
nlines=int(args.nl)
fit_type=str(args.type)
ident=str(args.i)
skipbot=totlines-nlines-skiptop
#Designates the file names
log_out=fit_type+'_fit.log'
file_out=file_d+'.'+fit_type+'_info.fit'
yfit=[]
perr=[]


#Defines functional forms
def linear(x, m, b):
    return m*x + b
def exp1(x, tau1, A):
    return A*np.exp(-x/tau1)
def exp2(x, tau1, tau2, A):
    return A*np.exp(-x/tau1) + (1-A)*np.exp(-x/tau2)
def exp3(x, tau1, tau2, tau3, A, B):
    return A*np.exp(-x/tau1) + B*np.exp(-x/tau2) + (1-A-B)*np.exp(-x/tau3)

#Reads the data file into xdata and ydata
xdata, ydata = np.genfromtxt(file_d, dtype=float, usecols=(col1,col2),autostrip=True,unpack=True,skip_header=skiptop, skip_footer=skipbot)



