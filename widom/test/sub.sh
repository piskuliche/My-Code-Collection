#MSUB -N lmps
#MSUB -q laird
#MSUB -d ./
#MSUB -j oe
#MSUB -l nodes=1:ppn=20:intel,mem=100gb,walltime=24:00:00

module purge
module load pylammps

mpirun -np 20 python widom_lammps.py in.ip 200 1000 > out
