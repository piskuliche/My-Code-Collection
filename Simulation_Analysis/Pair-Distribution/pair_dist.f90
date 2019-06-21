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
      integer :: ioerr, eweight

      real :: dr, L, rij_sq, pi
      real :: rho, const, r_hi, r_lo,ctmp
      real :: eavg
      real, dimension(3) :: rtmp, rij
      real, allocatable ::  h_id(:), tot_gofr(:), bl_gofr(:)
      real, allocatable :: tot_egofr(:), bl_egofr(:), en(:)
      real, allocatable :: r1(:,:,:), r2(:,:,:)
      real, allocatable :: g(:,:), eg(:,:)
      integer, allocatable :: typ(:), mol(:), mol1(:), mol2(:)
      integer, allocatable :: h(:,:)

      
      character(len=12) :: nfile, molfile, efile
      character(len=6) :: censtr, outstr,blstr
      
        
      ! Files Used:
      ! 11    Input: molfile
      ! 12    Input: nfile
      ! 13    Input: energies
      ! 20    Output: tot_gofr
      ! 21-41 Output: block gofr
      ! 42    Output: e_gofr
      ! 43-53 Output: block e_gofr


      pi = 3.14159
      
      NAMELIST /nml/ nfile, dr, startconfig, endconfig, startskip, endskip,&
      nblocks, selec1, selec2, L, or_int, natoms, molfile, eweight, efile

      ! Example Defaults
      nfile             = "traj.xyz"    ! Trajectory File Name
      molfile           = "molinfo.dat" ! Identifies molecules
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
      dr                = 0.1           ! Shell distance
      eweight           = 0             ! Weight by energies [0] off [1] on 
      efile             = "e_init.out"  ! Energy Weight File
      
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
      write(*,*) 'molfile: ', molfile
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
      write(*,*) 'dr: ', dr
      write(*,*) '*************'
    
      ! Set coords and bins
      dr = dr / L
      nb = FLOOR( 0.5/dr ) ! half box units
      nconfigs = (endconfig-startconfig)
      rconfigs = nconfigs/or_int
      write(*,*) "Starting Pair Distribution Calculation"
      write(*,*) "This will average over ", rconfigs, " configurations"
      write(*,*) "Opening ", nfile
      n1 = 0; n2 = 0
      ! Reads in the molecule info
      ! Sets the number of selec1 as n1
      ! Sets the number of selec2 as n2
      allocate(typ(natoms),mol(natoms))
      mol = 0; typ = 0
      open(11, file=trim(molfile), status='old')
      do j=1, natoms
        read(11,*) ctmp, mol(j), typ(j), ctmp
        if (typ(j) == selec1) then
            n1 = n1 + 1
        else if (typ(j) == selec2) then
            n2 = n2 + 1
        end if
      enddo
      if (selec1 == selec2) then
          n2 = n1
      end if
      close(11)
      write(*,*) "There are ",n1, " of selection 1"
      write(*,*) "There are ",n2, " of selection 2" 

      ! Allocate GofRArrays
      allocate(mol1(n1),mol2(n2)) ! Molecular Storage
      allocate(h_id(nb))          ! Ideal Result Storage
      allocate(h(rconfigs,nb), g(rconfigs,nb)) ! h and g arrays
      allocate(r1(rconfigs,n1,3), r2(rconfigs,n2,3)) ! coord arrays
      allocate(tot_gofr(nb),bl_gofr(nb)) ! gofr arrays
      ! Allocate Energy Fluctuations Arrays
      if (eweight == 1) then
        allocate(en(rconfigs))
        allocate(eg(rconfigs,nb)) ! eg array
        allocate(tot_egofr(nb), bl_egofr(nb)) ! Energy Flucts GofR
        en = 0.0; eg = 0.0; tot_egofr=0.0; bl_egofr=0.0
      endif

      mol1 = 0; mol2 = 0


      ! Calculations of ideal
      r_lo = 0.0; r_hi = 0.0; h_id = 0.0
      rho = REAL(n1) ! note that we are in box coords
      if (selec1 /= selec2) rho=REAL(n2)
      const = 4.0 * pi * rho / 3.0 ! Const for multiplying 
      do b=1, nb
        r_lo = REAL( b-1 ) * dr
        r_hi = r_lo + dr
        h_id(b) = const * ( r_hi**3.0 - r_lo ** 3.0 ) ! Ideal Number
      enddo

      ! Initialize to Zero
      h = 0; r1 = 0.0; r2 = 0.0; g = 0.0
      rij = 0.0; rij_sq = 0.0; tot_gofr=0.0

      ! Reads the trajectory file (and energy file)
      open(12, file=trim(nfile), status='old')
      if (eweight == 1) open(13, file=trim(efile), status='old')

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
        if (MOD(i,1000) == 0) then
            write(*,*) 'Configuration ',i,' read in.'
        end if
        ! Only operates if (i-1)%or_int is zero (allows skipping configs)
        if (MOD(i-1,or_int) == 0) then
            cnt = cnt + 1
            cnt1 = 1; cnt2 = 1
            if (eweight == 1) read(13,*) en(cnt) ! Reads in Energies
            do j=1, startskip
                read(12,*)
            enddo
            do j=1,natoms
                read(12,*) ctmp, (rtmp(a), a=1,3)
                ! Creates arrays of coordinates (in box coords)
                if (typ(j) == selec1) then
                    r1(cnt,cnt1,:) = rtmp(:)/L
                    mol1(cnt1) = mol(j)
                    ! Checks if selec1 and selec2 are the same
                    if(selec2 == selec1) then
                        r2(cnt,cnt1,:) = rtmp(:)/L
                        mol2(cnt1) = mol(j)
                    endif
                    cnt1 = cnt1 + 1
                ! Note - this will never happen if selec1 and selec2 are
                ! equivalent, as that is handeled previously.
                else if (typ(j) == selec2) then
                    r2(cnt,cnt2,:) = rtmp(:)/L
                    mol2(cnt2) = mol(j)
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
      if (eweight == 1) close(13)
      if (eweight == 1) eavg = SUM(en)/rconfigs
      write(*,*) "Calculating GofR"
      ! Loop over configurations
      !$OMP PARALLEL DO schedule(static) DEFAULT(NONE) &
      !$OMP& PRIVATE(i, j, k, rij, rij_sq, b, nb, n1, n2, dr, cnt) &
      !$OMP& SHARED(nconfigs, selec1, selec2, r1, r2, g, h, h_id,rconfigs) 
      do i=1, rconfigs
        do j=1, n1
            if (selec1 /= selec2) then
                do k=1, n2
                    if (mol1(j) /= mol2(k)) then
                        rij(:)      = r1(i,j,:) - r2(i,k,:)
                        rij(:)      = rij(:) - ANINT( rij(:) )
                        rij_sq      = SUM( rij**2 )
                        b           = FLOOR( SQRT( rij_sq ) / dr ) + 1
                        if ( b <= nb ) h(i,b) = h(i,b) + 1
                    endif
                enddo
            else
                if (j /= n1) then
                    do k=j+1, n1
                        if (mol1(j) /= mol2(k)) then
                            rij(:)  = r1(i,j,:) - r2(i,k,:)
                            rij(:)  = rij(:) - ANINT( rij(:) )
                            rij_sq  = SUM( rij**2 )
                            b       = FLOOR( SQRT( rij_sq ) / dr ) + 1
                            IF ( b <= nb ) h(i,b) = h(i,b) + 2
                        endif
                    enddo
                end if
            end if
        enddo
        ! Final Calculations
        g(i,:) = REAL( h(i,:) ) / REAL( n1 ) ! Norm by number of n1
        g(i,:) = g(i,:) / h_id(:) ! Normalize based on ideal result
        if (eweight == 1) eg(i,:) = en(i) * g(i,:) ! Weight by energies
      enddo
      !$OMP END PARALLEL DO
      write(*,*) "Tidying things up!"
      write(censtr, "(I0)") selec1
      write(outstr, "(I0)") selec2
      ! Total calculation
      do i=1, rconfigs
        do b=1,nb
            tot_gofr(b) = tot_gofr(b) + g(i,b)
            if (eweight == 1) tot_egofr(b) = tot_egofr(b) + eg(i,b)
        enddo
      enddo
      tot_gofr(:) = tot_gofr(:)/rconfigs
      if (eweight == 1) tot_egofr(:) = tot_egofr(:)/rconfigs - eavg*tot_gofr(:)
      ! Block calculation
      write(*,*) "Jiggying up blocks!"
      do j=1, nblocks
        bl_gofr=0.0
        do i=1,rconfigs/nblocks
            k=(j-1)*rconfigs/nblocks+i
            do b=1,nb
                bl_gofr(b) = bl_gofr(b) + g(k,b)
                if (eweight == 1) bl_egofr(b) = bl_egofr(b) + eg(k,b)
            enddo
        enddo
        bl_gofr(:)=bl_gofr(:)/(rconfigs/nblocks)
        if (eweight == 1) bl_egofr(:)=bl_egofr(:)/(rconfigs/nblocks)
        write(blstr, "(I0)") j
        open(20+j,file='bl_'//trim(blstr)//'_pairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
        if (eweight == 1) then
            open(42+j,file='bl_'//trim(blstr)//'_epairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
        endif
        do b=1,nb
            write(20+j,'(2f15.8)' ) (REAL(b)-0.5)*dr, bl_gofr(b)
            if (eweight == 1) then
                write(42+j,'(2f15.8)' ) (REAL(b)-0.5)*dr, bl_egofr(b)
            endif
        enddo
        close(20+j)
        if (eweight == 1) close(42+j)
      enddo


      ! Convert out of box units
      dr = dr*L
      
      write(*,*) "Converted back to original coordinates"
      ! Output to file
      write( *,* ) 'Outputing to output file'
      open(20, file='pairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
      if (eweight == 1) open(42,file='epairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
      do b=1, nb
        write(20,'(1f8.4,5x,1e15.6)' ) (REAL(b)-0.5)*dr, tot_gofr(b)
        if (eweight == 1) write(42,'(1f8.4,5x,1e15.6)' ) (REAL(b)-0.5)*dr, tot_egofr(b)
      enddo
      close(20)
      if (eweight == 1) close(42)
      write(*,*) "GofR Calculation Complete"
      if (eweight == 1) write(*,*) "Weighted GofR Calculation Complete"
     


END Program pairdist


