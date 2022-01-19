#!/usr/bin/env python

"""
Copyright 2022, Ezekiel Ashe Piskulich, Boston University

This is a small code that builds an empty input file and submission script for the grem_sort and grem_analyze fortran codes.

The only required input is the number of windows, passed as a command line argument.

"""

lipids={'DXPC':16,'DBPC':14,'DPPC':12,'DLPC':10,'DIPC':12,'DFPC':12,'DOPC':12,'POPC':12}
tails={10:"5 6 6 7 8 9 9 10", 12:"5 6 6 7 7 8 9 10 10 11 11 12", 14:"5 6 6 7 7 8 8 9 10 11 11 12 12 13 13 14", 16:"5 6 6 7 7 8 8 9 9 10 11 12 12 13 13 14 14 15 15 16"}

import sys
import numpy as np

nw = int(sys.argv[1])
nlip= int(sys.argv[2])
liptype=str(sys.argv[3])


f=open("grem.input",'w')

print("Creating grem.input")
f.write("#Nwindows\n")
f.write("%d\n" % nw)
f.write("#fstart fstop\n")
f.write("11 21\n")
f.write("#nlog dumpfreq\n")
f.write("10000 10\n")
f.write("#natoms nlipids\n")
f.write("%d %d\n" % (lipids[liptype],nlip))
f.write("#hatom\n")
f.write("2 2\n")
f.write("# lt1 lt2 lt3 lt4\n")
f.write("%s\n" % tails[lipids[liptype]]) 
print("Done")
f.close()

h=open("walkers.dat",'w')
for i in range(nw):
    h.write("walker-%d\n" % i)
h.close()

g=open("analyze.sh",'w')
g.write("#!/bin/bash\n#\n#$ -l h_rt=11:45:00\n#$ -j y\n#$ -N analyze\n#$ -pe omp 1\n#$ -V\n\n")

g.write("module load codecol/grem\n")
g.write("mkdir -p all_dumps\n")
g.write("grem_sort.exe\n\n\n")

g.write("cd all_dumps\n")
g.write("mv ../walkers.dat .\n")
g.write("for i in {0..%d};\n" % nw)
g.write("do\n")
g.write("    grem_analyze.exe <<< $i\n")
g.write("done\n")
g.close()




