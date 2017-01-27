#!/usr/bin/python

import math, sys
import numpy as np
import argparse
from argparse import RawTextHelpFormatter,SUPPRESS
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description='''Calculates the linear fit of data''',formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', help="File Name")
parser.add_argument('-c1', help='Numeric index of Column 0,1,2...N')
parser.add_argument('-c2', help='Numeric index of Column 0,1,2...N')
parser.add_argument('-s', help='# Lines of Data to skip')
parser.add_argument('-i',help='extra identifier')
args = parser.parse_args()
file_d=args.f
col1=int(args.c1)
col2=int(args.c2)
skip=int(args.s)
ident=str(args.i)

file_out=file_d+'.fit_info.fit'

log_out='linear_fit.log'
yfit = []
perr = []
def linear(x, m, b):
	return m*x + b
	

xdata,ydata = np.genfromtxt(file_d,dtype=float,usecols=(col1,col2),autostrip=True,unpack=True, skip_header=skip)
popt, pcov = curve_fit(linear, xdata, ydata, p0=(1,0))

yfit = linear(xdata,popt[0],popt[1])

n=len(xdata)
perr = np.sqrt(np.diag(pcov))
sum = open(str(ident+'lin_fit_sum.log'),'a')
sum.write('{0} {1} {2} {3}\n'.format(file_d,popt[0],popt[1],perr[0]))
log = open(log_out, 'a')
log.write('{0}\n'.format(file_d))
log.write('# Number of data points = {0}\n'.format(n))
log.write('# m = {0} b = {1}\n'.format(popt[0],popt[1]))
log.write('# m err = {0} b err = {1} \n'.format(perr[0],perr[1]))
log.write('# Std Dev = {0}\n'.format(np.sqrt(np.sum((linear(xdata,popt[0],popt[1])-ydata)**2)/n)))
np.savetxt(file_out, np.transpose([xdata,yfit]),fmt=['%10.5f','%10.5f'])
