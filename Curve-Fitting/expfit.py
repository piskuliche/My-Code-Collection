#!/usr/bin/python

import math, sys
import numpy as np
import argparse
from argparse import RawTextHelpFormatter,SUPPRESS
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description='''Calculates the exponential fit of data''',formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', help="File Name")
parser.add_argument('-c1', help='Numeric index of Column 0,1,2...N')
parser.add_argument('-c2', help='Numeric index of Column 0,1,2...N')
parser.add_argument('-s', help='# Lines of Data to skip')
parser.add_argument('-n', help='# number of exponentials to use')
parser.add_argument('-i', help='extra identifier')
args = parser.parse_args()
file_d=args.f
col1=int(args.c1)
col2=int(args.c2)
skip=int(args.s)
nexp=int(args.n)
ident=str(args.i)


log_out='exp_fit.log'
yfit = []
perr = []
def one(x, tau1, A):
    return A*np.exp(-x/tau1)
def two(x, tau1, tau2, A):
    return A*np.exp(-x/tau1) + (1-A)*np.exp(-x/tau2)
def three(x, tau1, tau2, tau3, A, B):
    return A*np.exp(-x/tau1) + B*np.exp(-x/tau2) + (1-A-B)*np.exp(-x/tau3)

xdata,ydata = np.genfromtxt(file_d,dtype=float,usecols=(col1,col2),autostrip=True,unpack=True, skip_header=skip)

#Single Exponential Fit
if nexp==1:
    #calls the single exponential function, fits it.
    popt, pcov = curve_fit(one, xdata, ydata, p0=(300,1))
    #defines the equation
    yfit = one(xdata,popt[0],popt[1])
    #defines the length of the data set
    n=len(xdata)
    #determines the standard deviation (aka sqrt(cov))
    perr = np.sqrt(np.diag(pcov))
    #opens the summary file
    sum = open(str(ident+'one_exp_smry.log'),'a')
    #writes the summary file (filename, tau1, A)
    sum.write('{0} {1} {2}\n'.format(file_d,popt[0],popt[1]))
    #opens the log file
    log = open(log_out, 'a')
    #writes the log file
    log.write('{0}\n'.format(file_d)) #filename
    log.write('# Number of data points = {0}\n'.format(n)) #number points
    log.write('# tau1 = {0} A = {1}\n'.format(popt[0],popt[1])) #fit parameters
    log.write('# tau1 err = {0} A err = {1} \n'.format(perr[0],perr[1])) #fit accuracy
    log.write('# Std Dev = {0}\n'.format(np.sqrt(np.sum((one(xdata,popt[0],popt[1])-ydata)**2)/n))) #STD DEV of fit
    file_out=file_d+".single.fitinfo"
    np.savetxt(file_out, np.transpose([xdata,yfit]),fmt=['%10.5f','%10.5f'])

#Double Exponential Fit
elif nexp==2:
    #calls the double exponential function, fits it.
    popt, pcov = curve_fit(two, xdata, ydata, p0=(5,40,.9))
    #defines the equation
    yfit = two(xdata,popt[0],popt[1],popt[2])
    #defines the length of the data set
    n=len(xdata)
    #determines the standard deviation (aka sqrt(cov))
    perr = np.sqrt(np.diag(pcov))
    #opens the summary file
    sum = open(str(ident+'two_exp_smry.log'),'a')
    #writes the summary file (filename, tau1, tau2, A)
    sum.write('{0} {1} {2} {3}\n'.format(file_d,popt[0],popt[1],popt[2]))
    #opens the log file
    log = open(log_out, 'a')
    #writes the log file
    log.write('{0}\n'.format(file_d)) #filename
    log.write('# Number of data points = {0}\n'.format(n)) #number points
    log.write('# tau1 = {0} tau2 = {1} A = {2}\n'.format(popt[0],popt[1],popt[2])) #fit parameters
    log.write('# tau1 err = {0} tau2 err = {1} A err = {2} \n'.format(perr[0],perr[1], perr[2])) #fit accuracy
    log.write('# Std Dev = {0}\n'.format(np.sqrt(np.sum((two(xdata,popt[0],popt[1],popt[2])-ydata)**2)/n))) #STD DEV of fit
    file_out=file_d+".double.fitinfo"
    np.savetxt(file_out, np.transpose([xdata,yfit]),fmt=['%10.5f','%10.5f'])
#Triple Exponential Function
else:
    #calls the triple exponential function, fits it.
    popt, pcov = curve_fit(three, xdata, ydata, p0=(50,400,1400,.2,.5))
    #defines the equation
    yfit = three(xdata,popt[0],popt[1],popt[2],popt[3],popt[4])
    #defines the length of the data set
    n=len(xdata)
    #determines the standard deviation (aka sqrt(cov))
    perr = np.sqrt(np.diag(pcov))
    #opens the summary file
    sum = open(str(ident+'three_exp_smry.log'),'a')
    #writes the summary file (filename, tau1, tau2, A)
    fullp=sorted([popt[0],popt[1],popt[2]])
    fullab=([popt[3],popt[4]])
    sum.write('{0} {1} {2} {3} {4} {5}\n'.format(file_d,fullp[0],fullp[1],fullp[2],fullab[0],fullab[1]))
    #opens the log file
    log = open(log_out, 'a')
    #writes the log file
    log.write('{0}\n'.format(file_d)) #filename
    log.write('# Number of data points = {0}\n'.format(n)) #number points
    log.write('# tau1 = {0} tau2 = {1} tau3 = {2} A = {3} B = {4}\n'.format(popt[0],popt[1],popt[2],popt[3],popt[4])) #fit parameters
    log.write('# tau1 err = {0} tau2 err = {1} tau3 err = {2} A err = {3} B err = {4}\n'.format(perr[0],perr[1], perr[2], perr[3], perr[4])) #fit accuracy
    log.write('# Std Dev = {0}\n'.format(np.sqrt(np.sum((three(xdata,popt[0],popt[1],popt[2],popt[3],popt[4])-ydata)**2)/n))) #STD DEV of fit
    #Single Exponential Fit
    file_out=file_d+".triple.fitinfo"
    np.savetxt(file_out, np.transpose([xdata,yfit]),fmt=['%10.5f','%10.5f'])

