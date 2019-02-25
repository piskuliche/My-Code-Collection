This is the diffusion calculation code. The fortran code can be run by the command

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
