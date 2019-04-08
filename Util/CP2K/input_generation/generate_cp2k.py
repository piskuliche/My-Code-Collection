####################################################################################################
# This is a code to create the genpot section of the CP2K input file
# It needs a simple file called param.in that includes the number of species and parameters.
#i.e.  the file should look like:
#---
#   4    nspecies
#   O   -0.8476   0.1553     3.166
#   H    0.4238   0.0000     0.000
#   Li   1.0000   0.3367     1.582
#   F   -1.0000   0.0074     4.514
#---
#If you do this - then you get a genpot section that you can pipe to an output file of your desire.
# Copyright Zeke Piskulich, University of Kansas, July 2018.
# You can contact me at piskuliche@ku.edu
####################################################################################################
import numpy as np

class param_init:
    def __init__(self, nsp):
        self.type = np.chararray((nsp,nsp),itemsize=2)
        self.type[:] = '  '
        self.q = np.zeros((nsp,nsp))
        self.eps = np.zeros((nsp,nsp))
        self.rm = np.zeros((nsp,nsp))

    def add_param(self, type, q, eps, rm, i):
        self.type[i,i] = type
        self.q[i,i] = q
        self.eps[i,i] = eps
        self.rm[i,i] = rm

    def mix_param(self,nsp):
        for i in range(nsp):
            for j in range(nsp):
                self.eps[i,j]=np.sqrt(self.eps[i,i]*self.eps[j,j])
                self.rm[i,j]=0.5*(self.rm[i,i]+self.rm[j,j])
    
    def print_genpot(self, nsp, rcut):
        for i in range(nsp):
            j=nsp-1
            while j >= i:
                print("\t&GENPOT")
                print("\t ATOMS %s %s" % (self.type[i,i],self.type[j,j]))
                print("\t FUNCTION D0*((R0/R)^12-2.0*(R0/R)^6")
                print("\t PARAMETERS D0 R0")
                print("\t UNITS kcalmol angstrom")
                print("\t VARIABLES R")
                print("\t VALUES %s %s" % (self.eps[i,j], self.rm[i,j]))
                print("\t RCUT %s" % rcut)
                print("\t&END GENPOT")
                j-=1

def generate():
    rcut = 12.0
    linecount = 0
    nspec=0
    with open('param.in', 'r') as inp:
        for line in inp:
            if linecount == 0:
                nspec = int(line.split()[0])
            linecount += 1
    param = param_init(nspec)
    linecount = 0
    with open('param.in', 'r') as inp:
        for line in inp:
            if linecount != 0:
               param.add_param(line.split()[0], line.split()[1], line.split()[2], line.split()[3],linecount-1)
            linecount += 1
    param.mix_param(nspec)
    param.print_genpot(nspec,rcut)
    
    
generate()
