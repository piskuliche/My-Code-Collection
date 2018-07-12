#This is an example code for the calculation of a Potential of Mean Force 
#from CP2K simulations.
#This sets up the various windows and gets them ready to run.
#Once you run this code - you should submit run_windows.sh
#Then you should run Alan Grossfield's WHAM code.
#That can be downloaded at http://membrane.urmc.rochester.edu/content/wham
#Copyright July 2018 Zeke Piskulich, University of Kansas.



import os, shutil
import numpy as np


start=1
end = 7
sep = end - start
bins = 20

meta = open('wham_metadata.info', 'w')
sub = open('run_windows.sh','w')

sub.write('#MSUB -N wham_bin\n')
sub.write('#MSUB -q sixhour\n')
sub.write('#MSUB -j oe\n')
sub.write('#MSUB -d ./\n')
sub.write('#MSUB -l nodes=1:ppn=10:intel,mem=100gb,walltime=6:00:00\n')
sub.write('#MSUB -t 0-%s\n\n\n\n' % (bins))

sub.write('module load cp2k/6.0/popt\n\n')


sub.write('cd $MOAB_JOBARRAYINDEX\n')
sub.write('mpirun -np 10 cp2k.popt inp_const.cp2k > run1.new\n')
sub.write("sed 's/\([ \t]\+[^ \t]*\)\{3\}$//' out.colvar > lif.distance\n")
sub.write("sed -i '1,20000 s/^/#/' lif.distance\n")
sub.write('cd ../\n')

for i in range(0,bins):
    ro = start+i*sep/float(bins)
    if not os.path.isdir(str(i)):
        os.mkdir(str(i))
    f = open(str(i)+'/collective.inc','w')
    f.write('\t&COLLECTIVE\n')
    f.write('\t COLVAR 1\n')
    f.write('\t INTERMOLECULAR TRUE\n')
    f.write('\t TARGET [angstrom] %s\n' % ro)
    f.write('\t\t&RESTRAINT\n')
    f.write('\t\t K 0.005\n')
    f.write('\t\t&END RESTRAINT\n')
    f.write('\t&END COLLECTIVE\n')
    f.close()
    shutil.copyfile('inp_const.cp2k',str(i)+'/inp_const.cp2k')

    meta.write('%s/lif.distance %s 0.01\n' % (str(i), ro))
    
meta.close()
sub.close

