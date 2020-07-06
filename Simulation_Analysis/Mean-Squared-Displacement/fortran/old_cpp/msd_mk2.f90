Program msd_rot_calc

    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, &
        error_unit, iostat_end, iostat_eor
    implicit none
    integer :: b, i, j, k, a, it, orig,bl
    integer :: nt, or_int, ind
    integer :: natms, nmols, atms_per_mol
    integer :: startskip, endskip, startconfig, endconfig, nconfigs, ntos
    integer :: nblocks, shamblock, nperblock
    integer :: msdindex, c2index1, c2index2
    integer :: com, counti
    integer :: tospblock
    real    :: dt, eavg, volume
    real    :: dp_O1, dp_O2, dsq
    real    :: MASS, cm_tmpmsd
   
    integer    :: ioerr
    logical    :: needblock

    real, dimension(3) :: L
    real, dimension(10) :: M, tmpmsd
    real, dimension(20) :: bl_cm_tmpmsd
    real, dimension(20,10) :: bl_tmpmsd
    real, dimension(100000,1000,10,3) :: r
    real, dimension(100000,1000,3) :: r_cm
    real, dimension(1000,10,3) :: r_old, shift
    real, dimension(1000,3) :: r_cm_old, shift_cm
    real, dimension(1000,3) :: e1, e2, e1_zero, e2_zero
    real, dimension(10000,0:5000,10) :: msd, c1, c2
    real, dimension(10000,0:5000) :: msd_cm
    

    character(len=10) nfile, mol_name
    character(len=2) ctmp
    character(len=2) blstr, shamstr

    NAMELIST /nml/ nfile, nt, or_int, dt, mol_name,nblocks &
        startskip, endskip, startconfig, endconfig, L, nmols, &
        needblock

    ! Example defaults
    nfile           = "traj.xyz" ! Traj File Name
    nt              = 400        ! How long the corr func is in dumps
    or_int          = 10         ! Separation of origins
    dt              = 0.01       ! Dump freq in ps
    L(1)            = 34.45      ! LX in Angstroms
    L(2)            = 34.45      ! LY in Angstroms
    L(3)            = 34.45      ! LZ in Angstroms
    mol_name        = "water"    ! molecule name
    nmols           = 400        ! nmols
    startconfig     = 0          ! Starting configuration
    endconfig       = 1000       ! Ending configuration
    startskip       = 2          ! Lines to skip at beginning of frame
    endskip         = 0          ! Lines to skip at end of frame
    nblocks         = 5          ! Default number of blocks
    needblock       = .true.     ! .true. turns on block averaging
    shamblock       = 6          ! acts as the psuedoblock for block aving
    
    ! Namelist from standard input
    READ ( unit=input_unit, nml=nml, iostat=ioerr)
    IF ( ioerr /= 0 ) THEN
        WRITE( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from stdin', ioerr
        IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
        IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of File'
        STOP 'Error in diffusion calc'
    END IF

    write( unit=output_unit, fmt='(a,t40,a)' ) ' File is ', nfile
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' Corr time is ', nt
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' Origin interval is ', or_int
    write( unit=output_unit, fmt='(a,t40,d)' ) ' Dump Freq is ', dt
    write( unit=output_unit, fmt='(a,t40,d)' ) ' Side length is ', L(1)
    write( unit=output_unit, fmt='(a,t40,a)' ) ' Mol Name is ', mol_name
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' startconfig is ', startconfig
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' endconfig is ', endconfig
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' startskip is ', startskip
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' endskip is ', endskip
    write( unit=output_unit, fmt='(a,t40,i15)' ) ' nblocks is ', nblocks
    write( unit=output_unit, fmt='(a,t40,L)' ) ' blockav on is ', needblock


    open(11,file=trim(mol_name)//'.txt',status='old')
    read(11,*) atms_per_mol
    read(11,*) msdindex
    read(11,*) c2index1, c2index2
    do i = 1, atms_per_mol
        read(11,*) M(i)
        MASS = MASS + M(i)
    enddo
    close(11)
    write( unit=output_unit,*) atms_per_mol 

    volume = L(1) * L(2) * L(3)
    nconfigs = endconfig-startconfig
    nperblock = nconfigs/nblocks

    ! Read the trajectory file
    write(unit=output_unit, fmt='(a,t40,i15)' ) ' startconfig is ', startconfig
    write(unit=output_unit, fmt='(a,t40,i15)' ) ' endconfig is ', endconfig

    ! Read the trajectory file
    open(12,file=trim(nfile),status='old') ! Open trajecotry file
    write( unit=output_unit,*) 'Reading file ', trim(nfile)
    
    ! skip to startconfig and throw out results
    do i=1, startconfig
        do j=1,startskip
            read(12,*)
        enddo
        do j=1, nmols
            do k=1, atms_per_mol
                read(12,*)
            enddo
        enddo
    enddo
    ! Loop over the blocks
    do b=1, nblocks
        write(*,*) 'Block ',i,' calculation'
        endconfig = startconfig + nperblock
        counti=1
        ! Read in the configurations
        do i=1, nperblock
            do j=1, startskip
                read(12,*)
            enddo
            do j=1, nmols
                do k=1, atms_per_mol
                    read(12,*) ctmp, (r(counti,a), a=1,3)
                    if ( k /= 1 ) then
                        ! Read in w/ PBC
                        r(counti,:) = r(counti,:)-L(:)*anint(( r(counti,:) - &
                            r(counti-k+1,:) )/L(:))
                    end if
                    ! Center of Mass Calc
                    r_cm(counti,:) += r(counti,:)*M(k)/MASS
                    counti += 1
                enddo
            enddo
            do j=1, endskip
                read(12,*)
            enddo
        enddo
        
        ! Writes config info to screen
        write(*,*) 'Block traj read complete.\n'
        write(*,*) 'There are ', nperblock, ' configurations in this block'

        ! Loop over time origins
        ind = 0
        do i=1, nconfigs-nt, or_int
            ! Loop over timesteps
            do it=i+2, nt+i
                do j-1, nmols
                    do k-1, atms_per_mol
                    ! Calculate the shift if it occurs
                    shift(:) = shift(:) - L(:)*anint((r(ind,:) - &
                        r_old(counti,:) )/L(:)
                    ! Calculate the square displacements
                    dsq = 
        startconfig += nperblock
    enddo
    close(12)
        


END Program msd_rot_calc 
