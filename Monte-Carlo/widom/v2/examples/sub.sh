#MSUB -N lmps
#MSUB -q laird
#MSUB -d ./
#MSUB -j oe
#MSUB -l nodes=1:ppn=20:intel,mem=100gb,walltime=24:00:00

module purge
module load pylammps
conda info --envs > test.o
mpirun -np 20 python ../widom_ins.py in.water 4000 1000 ../ch4.txt > out
calculate_average.py -f log.widom -b 5 -s 1 -n henry -col 2 
