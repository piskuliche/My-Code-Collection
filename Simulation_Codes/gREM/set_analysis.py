#!/usr/bin/env python
import sys
import numpy as np

nw = int(sys.argv[1])


f=open("grem.input",'w')

print("Creating grem.input")
f.write("#Nwindows\n")
f.write("%d\n" % nw)
f.write("#fstart fstop\n")
f.write("\n")
f.write("#nlog dumpfreq\n")
f.write("10000 10\n")
f.write("#natoms nlipids\n")
f.write("\n")
f.write("#hatom\n")
f.write("2 2\n")
f.write("# lt1 lt2 lt3 lt4\n")
print("Done")
f.close()

g=open("analyze.sh",'w')
g.write("#!/bin/bash\n#\n#$ -l h_rt=11:45:00\n#$ -j y\n#$ -N analyze\n#$ -pe omp 1\n#$ -V\n\n")

g.write("module load codecol/grem\n")
g.write("mkdir -p all_dumps\n")
g.write("grem_sort.exe\n\n\n")

g.write("cd all_dumps\n")
g.write("for i in {0..%d};\n" % nw)
g.write("do\n")
g.write("    grem_analyze.exe <<< $i\n")
g.write("done\n")
g.close()




