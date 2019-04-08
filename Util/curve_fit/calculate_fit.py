#!/usr/bin/python


#This code is a generalized fitting code which should be able to fit most space or tab delimited data. 
#
#Usage: python curve_fit.py -f filename -c1 0 -c1 1 -s1 0 -nl 100 -type exp3
#Copyright Zeke Piskulich, University of Kansas, 2017


import math, sys
import numpy as np
import argparse
from argparse import RawTextHelpFormatter,SUPPRESS
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
totlines=sum(1 for line in open(file_d))
skipbot=totlines-nlines-skiptop
#Designates the file names
log_out=fit_type+'_fit.log'
file_out=file_d+'.'+fit_type+'_info.fit'
summary=ident+fit_type+'_smry.log'
log=open(log_out, 'a')
smry=open(summary,'a')
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


#Where the fitting actually occurs
if fit_type == 'exp1' or fit_type =='EXP1':
    #calls the single exponential function
    popt, pcov = curve_fit(exp1, xdata, ydata, p0=(300,1))
    #defines the equations
    yfit = exp1(xdata, popt[0],popt[1])
    #determines the standard deviation
    perr = np.sqrt(np.diag(pcov))
    npar=2
    par=['tau1','A']
    stddev=np.sqrt(np.sum((exp1(xdata,popt[0],popt[1])-ydata)**2)/nlines)
elif fit_type =='exp2' or fit_type =='EXP2':
    #calls the double exponential function
    popt, pcov = curve_fit(exp2, xdata, ydata, p0=(5,2000,0.9))
    #defines the equations
    yfit = exp2(xdata, popt[0],popt[1], popt[2])
    #determines the standard deviation
    perr = np.sqrt(np.diag(pcov))
    npar=3
    par=['tau1','tau2','A']
    stddev=np.sqrt(np.sum((exp2(xdata,popt[0],popt[1],popt[2])-ydata)**2)/nlines)
elif fit_type =='exp3' or fit_type =='EXP3':
    #calls the triple exponential function
    popt, pcov = curve_fit(exp3, xdata, ydata, p0=(50,400,1400,0.2,0.5))
    #defines the equations
    yfit = exp3(xdata, popt[0],popt[1],popt[2],popt[3],popt[4])
    #determines the standard deviation
    perr = np.sqrt(np.diag(pcov))
    npar=5
    par=['tau1','tau2','tau3','A', 'B']
    stddev=np.sqrt(np.sum((exp3(xdata,popt[0],popt[1],popt[2],popt[3],popt[4])-ydata)**2)/nlines)
elif fit_type == 'linear' or fit_type == 'LINEAR':
    #calls the linear function
    popt, pcov = curve_fit(linear, xdata, ydata, p0=(1,0))
    #defines the equations
    yfit = linear(xdata, popt[0],popt[1])
    #determines the standard deviation
    perr = np.sqrt(np.diag(pcov))
    npar=2
    par=['m','b']
    stddev=np.sqrt(np.sum((linear(xdata,popt[0],popt[1])-ydata)**2)/nlines)
else:
    print "You have not selected a valid type of fit."
    raise SystemExit

#Writes the log file
log.write('{0}\n'.format(file_d))
log.write('# Number of Data Points: {0}'.format(nlines))
count=0
for item in par:
    log.write('# Parameter: {0} Error: {1}'.format(popt[count],perr[count]))
    count+=1
log.write('# STD DEV = {0}\n'.format(stddev))
np.savetxt(file_out, np.transpose([xdata,yfit]),fmt=['%10.5f','%10.5f'])

#Writes the smry file
smry.write('{0} '.format(file_d))
count=0
for item in par:
    smry.write('{0} '.format(popt[count]))
    count+=1
smry.write('\n')

