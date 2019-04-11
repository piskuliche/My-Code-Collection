#!/usr/bin/env python
"""This is a code that lets you farm many widom insertions to hopefully speed up the process"""
def write_slurm(inputfile, njobs):
    """Function that writes the sbatch script that submits the calculations"""
    f = open("run_widom.sh", 'w')

    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=farm_widom\n")
    f.write("#SBATCH --partition=laird,thompson,sixhour\n")
    f.write("#SBATCH --constraint=intel\n")
    f.write("#SBATCH --output=output_widom.log\n")
    f.write("#SBATCH --mail-user=piskuliche@ku.edu\n")
    f.write("#SBATCH --nodes=1\n")
    if njobs <= 20:
        f.write("#SBATCH --ntasks-per-node=%d\n" % njobs)
    else:
        f.write("#SBATCH --ntasks-per-node=20\n")
    f.write("#SBATCH --mem=%dG\n" % int(njobs/20.0*100))
    f.write("#SBATCH --time=6:00:00\n")
    f.write("module load codecol\n")
    f.write("module load lammps/2Dec2018\n")
    f.write("echo Time is `date`\n")
    f.write("echo Directory is `pwd`\n")
    f.write("\n")
    f.write('echo "Running on $SLURM_JOB_NODELIST nodes using $SLURM_CPUS_ON_NODE cores on each node"\n')
    f.write("\n")
    f.write("srun -n %d --multi-prog myrun.conf\n" % njobs)
    f.write("\n")
    f.write("echo Ending Time is `date`\n")
    f.write("exit 0\n")
    f.close()
    print("Slurm file written as run_widom.sh")
    return

def write_farms(inputfile, njobs):
    """Function that writes the farm script used in the submission"""
    f = open("myrun.conf", 'w')
    for i in range(njobs):
        f.write("%d widom_backbone.py -in %s -farm %s %%t\n" % (i,inputfile, njobs))
    f.close()
    print("Farm files written as myrun.conf")
    return

def write_lmps(njobs):
    """Function that writes rerun commands to new log files for each of the farmed trajectories"""
    for i in range(njobs):
        f = open("tmp_rerun"+str(i), 'w')
        f.write("log log.rerun%d\n" % i) 
        f.write("rerun dump_traj.widom%d dump x y z ix iy iz box yes\n" % i)
        f.close()
    return

def write_array(lmpsheader, njobs):
    f = open("run_ecalc.sh", 'w')

    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=ecalc_widom\n")
    f.write("#SBATCH --partition=laird,thompson,sixhour\n")
    f.write("#SBATCH --constraint=intel\n")
    f.write("#SBATCH --output=output_ecalc%A_%a.log\n")
    f.write("#SBATCH --mail-user=piskuliche@ku.edu\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --ntasks-per-node=4\n")
    f.write("#SBATCH --time=6:00:00\n")
    f.write("#SBATCH --mem=20G\n")
    f.write("#SBATCH --array=0-%d\n" % njobs)
    f.write("module load lammps/2Dec2018\n")
    f.write("cat %s tmp_rerun$SLURM_ARRAY_TASK_ID > in.rerun.widom$SLURM_ARRAY_TASK_ID\n" % (lmpsheader))
    f.write("mpirun lmp_mpi_opt -sf opt -in in.rerun.widom$SLURM_ARRAY_TASK_ID\n")
    f.close()




if __name__ == "__main__":
    import sys
    if "-h" in sys.argv or "-in" not in sys.argv or "-farm" not in sys.argv:
        print("This script creates a farming script for use with a SLURM Queue Engine")
        print("Usage python -in infile -farm numberjobs -header lmpsheader")
        sys.exit()

    if "-in" in sys.argv:
        index = sys.argv.index("-in")+1
        inputfile = str(sys.argv[index])
    else:
        print("Need to include an inputfile argument")
        sys.exit()

    if "-farm" in sys.argv:
        index = sys.argv.index("-farm")+1
        farmjobs = int(sys.argv[index])
    else:
        print("Need to include a farm argument")
        sys.exit()

    if "-header" in sys.argv:
        index = sys.argv.index("-header")+1
        head  = str(sys.argv[index])
    else:
        print("Need to include a header argument")
        sys.exit()

    print("Writing the necessary files")
    write_slurm(inputfile, farmjobs)
    write_farms(inputfile, farmjobs)
    write_lmps(farmjobs)
    write_array(head,farmjobs)
    print("Farm script complete")
    print("Steps")
    print("1) submit run_widom.sh")
    print("2) submit run_ecalc.sh")
    print("3) interactively run the following")
    print("widom_calculation.py -in %s -farm %d" % (inputfile, farmjobs))
    print("Have a good day!")

    

    
