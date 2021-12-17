#!/usr/bin/env python
import numpy as np
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-start', default=0,type=int, help='Starting index for calculating [default=0]')
parser.add_argument('-stop', default=5, type=int, help='Ending index for calculating [default=5]')
parser.add_argument('-nbins', default=100, type=int, help='Number of bins [default=100]')
parser.add_argument('-eta', default=-0.03, type=float, help='Eta value used in the simulations [default=-0.03]')
parser.add_argument('-Ho', default=-30000, type=float, help='H0 value used in the simulations [default=-30000]')
parser.add_argument('-nb', default=5, type=int, help='Number of blocks used for block averaging [default=5]')
parser.add_argument('-nprocs', default=4, type=int, help='Number of processors [default=4]')
parser.add_argument('-wprocs', default=16, type=int, help='Number of walkdown processors [default=16]')
parser.add_argument('-name', default="dppc", type=str, help='Name of molecule')
parser.add_argument('-hatom', default=4, type=int, help='Atom type for header atom')
parser.add_argument('-rcut', default=12, type=float, help='Cutoff distance in A')
args= parser.parse_args()


startval=args.start
stopval=args.stop
nbins = args.nbins
eta=args.eta
Ho=args.Ho
nblocks=args.nb
nprocs=args.nprocs
wprocs=args.wprocs
hatom=args.hatom
rcut=args.rcut
name=args.name

with open("autogen.log",'w') as log:
    for item in sys.argv:
        log.write("%s " % item)


with open("step_setup.sh",'w') as f:
    f.write("#!/bin/bash\n")
    f.write("module load codecol/grem\n")
    f.write("\n")
    # Run Setup
    f.write("cp ../build/equil/200/final_restart_file ./restart_file\n")
    f.write("cp ../build/equil/system.in.* .\n")
    f.write("cp -r ../../base .\n")
    f.write("echo dielectric 15.0 >> system.in.init\n")
    f.write("run_gREM.py -setup 1 -start %d -stop %d -nbins %d\n" % (startval, stopval, nbins))
    # Check restart_file is in main directory
    f.write("if [ ! -f 'restart_file' ]; then echo 'Error: restart_file does not exist';exit; fi\n")
    # Copy restart_file to each directory
    spacing=(stopval-startval)/(nbins-1)
    f.write("for i in {%d..%d}; do cp restart_file $i; done\n" % (0,nbins-1))

with open("step_walkdown.sh",'w') as f:
    f.write("#!/bin/bash\n")
    f.write("module load codecol/grem\n")
    f.write("\n")
    f.write("if [ ! -f 'base/walkdown.in' ]; then echo 'Error: base/walkdown.in does not exist';exit; fi\n")
    f.write("cp base/walkdown.in .\n")
    f.write("if [ ! -f 'grem.include' ]; then echo 'Error: grem.include does not exist';exit; fi\n")
    f.write("if [ ! -f 'system.in.settings' ]; then echo 'Error: system.in.settings does not exist';exit; fi\n")
    # Setup walkdown script
    f.write("run_gREM.py -start %d -stop %d -eta %5.2f -Ho %8.2f -nprocs %d -inp walkdown.in -lambdafile grem.include -walkdown 1\n" % (0,nbins,eta,Ho,wprocs))
    f.write('sed -i -e "s@dppc@%s@g" walkdown.sh\n' % name)
    f.write('qsub walkdown.sh\n')

with open("step_grem.sh",'w') as f:
    f.write("#!/bin/bash\n")
    f.write("module load codecol/grem\n")
    f.write("\n")
    f.write("cp base/start.gREM .\n")
    f.write("cp base/run.gREM-TEMPLATE . \n")
    f.write("cp base/submit-temper.bash .\n")
    # Update 
    f.write("sed -i -e 's@AAA@%5.2f@g' start.gREM\n" % eta)
    f.write("sed -i -e 's@AAA@%5.2f@g' run.gREM-TEMPLATE\n" % eta)
    f.write("sed -i -e 's@BBB@%8.2f@g' start.gREM\n" % Ho)
    f.write("sed -i -e 's@BBB@%8.2f@g' run.gREM-TEMPLATE\n" % Ho)
    f.write("sed -i -e 's@AAA@%s@g' submit-temper.bash\n" % name)
    usedprocs=nprocs*nbins
    num_node = np.ceil(usedprocs/16.0)
    num_cores = num_node*16
    f.write("sed -i -e 's@BBB@%d@g' submit-temper.bash\n" % num_cores)
    f.write("sed -i -e 's@CCC@%d@g' submit-temper.bash\n" % usedprocs)
    f.write("sed -i -e 's@DDD@%d@g' submit-temper.bash\n" % nbins)

with open("step_analyze.sh",'w') as f:
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("mkdir -p anal\n")
    f.write("cd anal\n")
    f.write("run_gREM.py -start 1 -stop 11 -eta %5.2f -Ho %8.2f -workdir .. -nbins 100 -nb %d -restart 0\n" % (eta, Ho, nblocks))
    f.write("run_gREM.py -start 1 -stop 11 -eta %5.2f -Ho %8.2f -workdir .. -nbins 100 -nb %d -restart 1\n" % (eta,Ho,nblocks))
    f.write("cd ..\n")

