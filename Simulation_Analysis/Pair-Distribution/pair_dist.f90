Program pairdist
      ! Things left to be done: 
        ! Need to bug check
      Use, INTRINSIC :: iso_fortran_env, ONLY: input_unit, output_unit, &
        error_unit, iostat_end, iostat_eor
      implicit none
      integer :: i, j, k, a, b             ! Loop iteration variables
      integer :: cnt, cnt1, cnt2, nb, tmp    ! Counters
      integer :: selec1, selec2, n1, n2 ! Type Selections, number each type
      integer :: nconfigs, startconfig, endconfig, rconfigs
      integer :: startskip, endskip
      integer :: natoms, or_int, nblocks
      integer :: ioerr

      real :: dr, L, rij_sq, pi
      real :: rho, const, r_hi, r_lo
      real, dimension(3) :: rtmp, rij
      real, allocatable :: typ(:), h_id(:), tot_gofr(:), bl_gofr(:)
      real, allocatable :: r1(:,:), r2(:,:)
      real, allocatable :: g(:,:)
      integer, allocatable :: h(:,:)

      
      
      character(len=10) :: nfile
      character(len=2) :: selstr1, selstr2
      
      pi = 3.14159
      
      NAMELIST /nml/ nfile, dr, startconfig, endconfig, startskip, endskip,&
      nblocks, selec1, selec2, L

      ! Example Defaults
      nfile             = "traj.xyz"    ! Trajectory File Name
      L                 = 34.0          ! Simulation box length (angstrom)
      natoms            = 1000          ! Number of Atoms
      startconfig       = 0             ! Starting Configuration (Frame)
      endconfig         = 1000          ! Ending Configuration (Frame)
      or_int            = 1000          ! Number of frames between config
      startskip         = 2             ! Number lines to skip, beg of frame
      endskip           = 0             ! Number lines to skip, end of frame
      nblocks           = 5             ! Number of blocks, block average
      selec1            = 1             ! Atom Type 1 for GofR
      selec2            = 1             ! Atom Type 2 for GofR
      
      ! Reads the namelist from standard input
      read ( unit=input_unit, nml=nml, iostat=ioerr)
      if ( ioerr /= 0 ) then
        write( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from stdin', ioerr
        if ( ioerr == iostat_eor ) write( unit=error_unit, fmt='(a)') 'End of record'
        if ( ioerr == iostat_end ) write( unit=error_unit, fmt='(a)') 'End of file'
      end if
      write(*,*) '*************'
      write(*,*) '*   START   *'
      write(*,*) '*************'
      write(*,*) 'Input Details'
      write(*,*) '*************'
      write(*,*) 'nfile: ', nfile
      write(*,*) 'L: ', L
      write(*,*) 'natoms: ', natoms
      write(*,*) 'startconfig: ', startconfig
      write(*,*) 'endconfig: ', endconfig
      write(*,*) 'or_int: ', or_int
      write(*,*) 'startskip: ', startskip
      write(*,*) 'endskip: ', endskip
      write(*,*) 'nblocks: ', nblocks
      write(*,*) 'selec1: ', selec1
      write(*,*) 'selec2: ', selec2
      write(*,*) '*************'
      
      ! Set coords and bins
      dr = dr / L
      nb = FLOOR( 0.5/dr ) ! half box units
      nconfigs = (endconfig-startconfig)
      rconfigs = nconfigs/or_int

      ! read 1 frame and calc n1, n2
      open(19, file=trim(nfile), status='old')
      do j=1, startskip
        read(19,*)
      enddo
      do j=1, natoms
        read(19,*) tmp, (rtmp(a), a=1,3)
        if (tmp == selec1) then
            n1 = n1 + 1
        else if (tmp == selec2) then
            n2 = n2 + 1
        end if
      enddo
      if (selec1 == selec2) then
          n2 = n1
      end if
      close(19)

      ! Allocate Arrays
      allocate(typ(natoms))
      allocate(h_id(nb))
      allocate(h(rconfigs,nb), g(rconfigs,nb))
      allocate(r1(rconfigs,n1), r2(rconfigs,n2))

      ! Calculations of ideal
      r_lo = 0.0; r_hi = 0.0; h_id = 0.0
      rho = REAL(n1) ! note that we are in box coords
      const = 4.0 * pi * rho / 3.0 ! Const for multiplying 
      do b=1, nb
        r_lo = REAL( b-1 ) * dr
        r_hi = r_lo + dr
        h_id(b) = const * ( r_hi**3.0 - r_lo ** 3.0 ) ! Ideal Number
      enddo

      ! Initialize to Zero
      h = 0; r1 = 0.0; r2 = 0.0; g = 0.0
      rij = 0.0; rij_sq = 0.0

      ! Reads the trajectory file 
      open(12, file=trim(nfile), status='old')
      write( unit=output_unit, *) 'Reading file ', trim(nfile)
      do i=1, startconfig
        do j=1, startskip
            read(12,*)
        enddo
        do j=1, natoms
            read(12,*)
        enddo
      enddo
      cnt = 0
      do i=1, nconfigs
        if (MOD(i,5000) == 0) then
            write(*,*) 'Configuration ',i,' read in.'
        end if
        ! Only operates if (i-1)%or_int is zero (allows skipping configs)
        if (MOD(i-1,or_int) == 0) then
            cnt = cnt + 1
            do j=1, startskip
                read(12,*)
            enddo
            do j=1,natoms
                read(12,*) typ(j), (rtmp(a), a=1,3)
                ! Creates arrays of coordinates (in box coords)
                if (typ(j) == selec1) then
                    r1(cnt1,:) = rtmp(:)/L
                    ! Checks if selec1 and selec2 are the same
                    if(selec2 == selec1) then
                        r2(cnt1,:) = rtmp(:)/L
                    endif
                    cnt1 = cnt1 + 1
                ! Note - this will never happen if selec1 and selec2 are
                ! equivalent, as that is handeled previously.
                else if (typ(j) == selec2) then
                    r2(cnt2,:) = rtmp(:)/L
                    cnt2 = cnt2 + 1
                end if
            enddo
        else
            do j=1, startskip
                read(12,*)
            enddo
            do j=1, natoms
                read(12,*)
            enddo
        end if
      enddo
      close(12)
      ! Loop over configurations
      !$OMP PARALLEL DO schedule(static) DEFAULT(NONE) &
      !$OMP& PRIVATE(i, j, k, rij, rij_sq, b, nb, n1, n2, dr, cnt) &
      !$OMP& SHARED(nconfigs, selec1, selec2, r1, r2, g, h, h_id) 
      do i=1, nconfigs
        do j=1, n1
            if (selec1 /= selec2) then
                do k=1, n2
                    rij(:)      = r1(j,:) - r2(k,:)
                    rij(:)      = rij(:) - ANINT( rij(:) )
                    rij_sq      = SUM( rij**2 )
                    b           = FLOOR( SQRT( rij_sq ) / dr ) + 1
                    IF ( b <= nb ) h(i,b) = h(i,b) +2
                enddo
            else
                if (i /= n1) then
                    do k=j+1, n1
                        rij(:)  = r1(j,:) - r2(k,:)
                        rij(:)  = rij(:) - ANINT( rij(:) )
                        rij_sq  = SUM( rij**2 )
                        b       = FLOOR( SQRT( rij_sq ) / dr ) + 1
                        IF ( b <= nb ) h(i,b) = h(i,b) + 2
                    enddo
                end if
            end if
        enddo
        ! Final Calculations
        g(i,:) = REAL( h(i,:) ) / REAL( n1 * cnt )
        g(i,:) = g(i,:) / h_id(:)
      enddo
      !$OMP END PARALLEL DO

      ! Block calculation
        

      ! Convert out of box units
      dr = dr*L

      ! Output to file
      write(selstr1,"(I0)") selec1
      write(selstr2,"(I0)") selec2
      write( unit=output_unit, * ) 'Outputing to output file'
      open(13, file='pairdist_'//trim(selstr1)//'_'//trim(selstr2)//'_'//'.dat')
      do b=1, nb
        write(13,'(2f15.8)' ) (REAL(b)-0.5)*dr, tot_gofr(b)
      enddo
      close(13)
     


END Program pairdist


