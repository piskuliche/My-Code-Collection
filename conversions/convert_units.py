#!/bin/bash/env python
import sys

class constants:
    def __init__(self):
        self.Na  = 6.02214E23       #1/mol
        self.F   = 96485.33         #C/mol
        self.amu = 1.660538E-27     # kg
        self.R   = 8.3144           # J/(mol K)
        self.ke  = 8.987551E9       # N m^2 / C^2
        self.c   = 299792458        # m/s
        self.kb  = 1.38065E-23      # J/K
        self.e   = 1.602176E-19     # C
        self.g   = 9.80665          # m/s^2
        self.Rin = 1.0973731568E7   # 1/m
        self.h   = 6.62607004E-34   # J s
    def printall(self):
        print('Na  %s 1/mol' % self.Na)
        print('F   %s C/mol' % self.F)
        print('amu %s kg' % self.amu)
        print('R   %s J/(mol K)' % self.R)
        print('ke  %s N m^2 / C^2' % self.ke)
        print('c   %s m/s' % self.c)
        print('kb  %s J/K' % self.kb)
        print('e   %s C' % self.e)
        print('g   %s m/s^2' % self.g)
        print('Rin %s 1/m' % self.Rin)
        print('h   %s J s' % self.h)

n = constants()

n.printall()
