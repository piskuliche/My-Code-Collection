Program msd_rot_calc

    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, &
        error_unit, iostat_end, iostat_eor
    implicit none
    integer :: i, j, k, a, it, orig,bl
    integer :: nt, or_int
    integer :: natms, nmols, atms_per_mol
    integer :: startskip, endskip, startconfig, endconfig, nconfigs, ntos
    integer :: nblocks
    integer :: msdindex, c2index1, c2index2
    integer :: com, counti
    integer :: tospblock
    real    :: dt, eavg, volume
    real    :: dp_O1, dp_O2, dsq
    real    :: MASS, cm_tmpmsd
   
    integer    :: ioerr

    real, dimension(3) :: L
    real, dimension(10) :: M, tmpmsd
    real, dimension(20) :: bl_cm_tmpmsd
    real, dimension(20,10) :: bl_tmpmsd
    real, dimension(500000,1000,10,3) :: r
    real, dimension(500000,1000,3) :: r_cm
    real, dimension(1000,10,3) :: r_old, shift
    real, dimension(1000,3) :: r_cm_old, shift_cm
    real, dimension(1000,3) :: e1, e2, e1_zero, e2_zero
    real, dimension(10000,0:5000,10) :: msd, c1, c2
    real, dimension(10000,0:5000) :: msd_cm
    

    character(len=10) nfile, mol_name
    character(len=2) ctmp

    NAMELIST /nml/ nfile, nt, or_int, dt, mol_name,nblocks &
        startskip, endskip, startconfig, endconfig, L, nmols

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
    
    
    ! Read the trajectory file
    open(12,file=trim(nfile),status='old') ! Open trajecotry file
    write( unit=output_unit,*) 'Reading file ', trim(nfile)
    
    do i=1, startconfig
        do j=1, startskip
            read(12,*)
        enddo
        do j=1, nmols
            do k=1, atms_per_mol
                read(12,*)
            enddo
        enddo
        do j=1,endskip
            read(12,*)
        enddo
    enddo
    nconfigs = endconfig-startconfig
    
    do i=1, endconfig-startconfig
        if (MOD(i,5000) == 0) then
            write(*,*) 'Configuration ',i,'read in.'
        end if
        do j=1, startskip
            read(12,*)
        enddo
        do j=1, nmols
            do k=1, atms_per_mol
                read(12,*) ctmp, (r(i,j,k,a), a=1,3)
                if ( k /= 1 ) then
                    r(i,j,k,:) = r(i,j,k,:) - L(:)*anint(( r(i,j,k,:) - &
                        r(i,j,1,:) )/L(:))
                end if
                r_cm(i,j,:) = r_cm(i,j,:)+r(i,j,k,:)*M(k)/MASS
            enddo
        enddo
        do j=1, endskip
            read(12,*)
        enddo
    enddo
    close(12)

    write(unit=output_unit,*) 'Trajectory Read Complete.\n'
    write(unit=output_unit,*) 'There are ', nconfigs, ' configurations.'

    ! Loop over time origins
    counti = 0
    do i=1, nconfigs-nt, or_int
        if(MOD(i*or_int,5000) == 0) then
            write(*,*) 'Reached the ', i,'th time origin'
        end if
        ! Set the Old Coordinates
        counti = counti + 1
        r_old(:,:,:) = r(i,:,:,:)
        r_cm_old(:,:) = r_cm(i,:,:)
        shift = 0.0
        shift_cm = 0.0

        ! Loop over the timesteps in each trajectory
        do it=i+2, nt+i
            ! Loop over molecules
            do j = 1, nmols
                do k=1, atms_per_mol
                    ! Calculate the shift if it occurs.
                    shift(j,k,:) = shift(j,k,:) - L(:)*anint((r(it,j,k,:) - &
                        r_old(j,k,:) )/L(:))
                    ! Calculate the square displacements
                    dsq = ( r(it,j,k,1) + shift(j,k,1) - r(i,j,k,1) ) ** 2. &
                         +( r(it,j,k,2) + shift(j,k,2) - r(i,j,k,2) ) ** 2. &
                         +( r(it,j,k,3) + shift(j,k,3) - r(i,j,k,3) ) ** 2.
                    msd(counti, it-1-i, k) = msd(counti, it-1-i, k) + dsq
                    ! Calculate the contribution to the c1,c2
                enddo ! End Atoms Loop (k)
                ! Calculate the shift if it occurs.
                shift_cm(j,:) = shift_cm(j,:) - L(:)*anint((r_cm(it,j,:) - &
                                r_cm_old(j,:) )/L(:))
                ! Calculate the square displacements
                dsq = ( r_cm(it,j,1) + shift_cm(j,1) - r_cm(i,j,1) ) ** 2. &
                    +( r_cm(it,j,2) + shift_cm(j,2) - r_cm(i,j,2) ) ** 2. &
                    +( r_cm(it,j,3) + shift_cm(j,3) - r_cm(i,j,3) ) ** 2.
                msd_cm(counti,it-1-i) = msd_cm(counti, it-1-i) + dsq
            enddo ! End Molecules Loop (j)
            r_old(:,:,:) = r(it,:,:,:)
            r_cm_old(:,:) = r_cm(it,:,:)
       enddo ! End t's loop (it)
    enddo ! End Time origin loop (i)
    
    ntos = (nconfigs-nt)/real(or_int)
    write(*,*) ntos, nconfigs, nt, or_int
    tmpmsd=0.0
    cm_tmpmsd=0.0
    bl_tmpmsd=0.0
    bl_cm_tmpmsd=0.0
    open(20,file='msd_'//trim(mol_name)//'.dat') ! open MSD output file
    open(21,file='cm_msd_'//trim(mol_name)//'.dat') ! open COM MSD output file
    do it=1,nt
        do i=1, ntos
            do k=1,atms_per_mol
                tmpmsd(k) = tmpmsd(k) + msd(i, it-1, k)
            enddo
        cm_tmpmsd = cm_tmpmsd + msd_cm(i, it-1)
        enddo
        do bl=1, nblocks
            tospblock = ntos/nblocks
            do i=bl*tospblock-tospblock+1, bl*tospblock+1
                do k=1,atms_per_mol
                    bl_tmpmsd(bl,k) = bl_tmpmsd(bl,k) + msd(i, it-1, k)
                enddo
                bl_cm_tmpmsd(bl) = bl_cm_tmpmsd(bl) + msd_cm(i,it-1)
            enddo 
        enddo
        write(20,'(100f12.5)') real(it*dt), (tmpmsd(k)/real(nmols)/real(ntos), (bl_tmpmsd(bl,k)/real(nmols)/real(ntos)*nblocks, bl=1,nblocks), k=1,atms_per_mol)
        write(21,'(100f12.5)') real(it*dt), cm_tmpmsd/real(nmols)/real(ntos), (bl_cm_tmpmsd(bl)/real(nmols)/real(ntos)*nblocks, bl=1,nblocks) 
        do k=1,10
            tmpmsd(k) = 0.0
            bl_tmpmsd(:,k)=0.0
        enddo
        cm_tmpmsd = 0.0
        bl_cm_tmpmsd(:)=0.0
    enddo
    close(20)
    close(21)
    

                


   



END Program msd_rot_calc 
