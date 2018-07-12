#MSUB -N cp2k-run
#MSUB -q laird
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=20:intel,mem=100gb,walltime=24:00:00

module load cp2k/6.0/popt

echo Time is `date`
echo Directory is `pwd`

mpirun -np 20 cp2k.popt inp_const.cp2k > run1.new

echo Ending Time is `date`
exit 0
