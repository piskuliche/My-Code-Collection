This is the diffusion calculation code. The fortran code can be run by the command

Compile with 
 module load compiler/pgi/18
 pgf90 -O3 main.f90 -mcmodel=medium -o calc_msd.exe

or: 
 module load compiler/pgi/19
 pgfortran -O3 main_pgi.f90 -mp -mcmodel=medium -o calc_msd.exe
 which compiles it in parallel.

calc_msd.exe < diff.inp

where diff.inp has the following form:

&nml
nt=400,
or_int=10,
dt=0.5,
/

! Example defaults
    nfile           = "traj.xyz" ! Traj File Name
    nt              = 400        ! How long the corr func is in dumps
    or_int          = 10         ! Separation of origins
    dt              = 0.01       ! Dump freq in ps
    L(1)            = 34.45      ! Volume in Angstroms
    L(2)            = 34.45      !
    L(3)            = 34.45
    mol_name        = "water"    ! molecule name
    nmols           = 400        ! nmols
    startconfig     = 0          ! Starting configuration
    endconfig       = 1000       ! Ending configuration
    startskip       = 2          ! Lines to skip at beginning of frame
    endskip         = 0          ! Lines to skip at end of frame
    nblocks         = 5          ! Default number of blocks
    needblock       = .true.     ! Defaults to running blocks

This code is enabled for openmp, if you want to use this with openmp do the following:
1) compile the code with the mp flag as pgfortran -O3 -mp main.f90 -mcmodel=medium -o calc_msd.exe
2) when running the code interactively you need to set two environment variables - I am updating the modulefile to set these automatically

    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK (or a number)
    export OMP_STACKSIZE=3G

3) To run these within a sbatch script include the two above settings and run as:

    srun  --unbuffered calc_msd.exe n < foo.inp

    where n is only included if needblock  is .false. and is equal to shamblock which goes from 0 to nblocks-1



