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
      integer :: ioerr, eweight, reweight

      real*8 :: dr, L, rij_sq, pi, minL
      real*8 :: rho, const, r_hi, r_lo,ctmp
      real*8 :: eavg, desqavg,decubavg, weight, sumweight
      real*8, dimension(3) :: rtmp, rij
      real*8, allocatable :: eavgbl(:), desqavgbl(:), decubavgbl(:)
      real*8, allocatable ::  h_id(:), tot_gofr(:), bl_gofr(:)
      real*8, allocatable :: tot_egofr(:), bl_egofr(:), tot_e2gofr(:), tot_e3gofr(:), bl_e2gofr(:),bl_e3gofr(:), en(:),den(:)
      real*8, allocatable :: r1(:,:,:), r2(:,:,:)
      real*8, allocatable :: g(:,:), eg(:,:), e2g(:,:), e3g(:,:), npth_id(:,:)
      real*8, allocatable :: nptL(:), nptdr(:), nptconst(:)
      integer, allocatable :: typ(:), mol(:), mol1(:), mol2(:)
      integer, allocatable :: h(:,:)
      logical :: file_exists
      
      character(len=20) :: nfile, molfile, efile
      character(len=6) :: censtr, outstr, blstr
      character(len=6) :: etype
      
        
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
      nblocks, selec1, selec2, L, or_int, natoms, molfile, eweight, efile, &
      etype, reweight

      ! Example Defaults
      nfile             = "traj.xyz"    ! Trajectory File Name
      molfile           = "molinfo.dat" ! Identifies molecules
      L                 = 30.0          ! Simulation box length (angstrom)
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
      etype             = "e"           !
      reweight          = 0             ! Should reweight by prob dist [0] off [1] on

      
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
      if ( L > 0.0) write(*,*) 'L: ', L
      if ( L == 0.0) write(*,*) 'L: NPT'
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


      ! Set number of configurations
      nconfigs = (endconfig-startconfig)
      rconfigs = nconfigs/or_int

      ! Read box length if NPT
      if ( L == 0.0 ) then
          inquire(file='L.dat', exist=file_exists)
          write(*,*) "L.dat exists: ", file_exists
          open(15, file='L.dat', status='old')
          allocate(nptL(rconfigs), nptdr(rconfigs), nptconst(rconfigs))
          nptL=0.0; cnt = 0; minL = 1000.0; nptdr=0.0; nptconst=0.0
          do i=1, nconfigs
              if (MOD(i-1,or_int) == 0) then
                  cnt = cnt + 1
                  read(15,*) nptL(cnt)
                  ! Sets minimum length
                  if (nptL(cnt) < minL) minL=nptL(cnt)
                  ! calculates the roving dr
                  nptdr(cnt) = dr/nptL(cnt)
              endif
              ! Calculates number of bins based on minimum half box size of run.
              nb = FLOOR( 0.5/(dr/minL) )
          enddo
      else
          dr = dr / L
          nb = FLOOR( 0.5/dr ) ! half box units
      end if
    
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
      allocate(h_id(nb),npth_id(rconfigs,nb))          ! Ideal Result Storage
      allocate(h(rconfigs,nb), g(rconfigs,nb)) ! h and g arrays
      allocate(r1(rconfigs,n1,3), r2(rconfigs,n2,3)) ! coord arrays
      allocate(tot_gofr(nb),bl_gofr(nb)) ! gofr arrays
      allocate(eavgbl(nblocks)) ! block energy array
      ! Allocate Energy Fluctuations Arrays
      if (eweight == 1) then
        allocate(en(rconfigs),den(rconfigs))
        allocate(eg(rconfigs,nb)) ! eg array
        allocate(e2g(rconfigs,nb)) !e2g array
        allocate(e3g(rconfigs,nb)) !e3g array
        allocate(tot_egofr(nb), bl_egofr(nb)) ! Energy Flucts GofR
        allocate(tot_e2gofr(nb), bl_e2gofr(nb)) ! E2 Flucts GofR
        allocate(tot_e3gofr(nb), bl_e3gofr(nb)) ! E3 Flucts GofR)
        allocate(desqavgbl(nblocks), decubavgbl(nblocks))
        en = 0.0; den = 0.0; eg = 0.0; tot_egofr=0.0; bl_egofr=0.0
        tot_e2gofr=0.0; bl_e2gofr=0.0; desqavg=0.0; e2g=0.0
        tot_e3gofr=0.0; bl_e3gofr=0.0; decubavg=0.0; e3g=0.0
        
      endif

      mol1 = 0; mol2 = 0


      ! Calculations of ideal
      r_lo = 0.0; r_hi = 0.0; h_id = 0.0; npth_id =0.0
      rho = REAL(n1) ! note that we are in box coords
      if (selec1 /= selec2) rho=REAL(n2)

      ! Sets the constants for an nvt run
      const = 4.0 * pi * rho / 3.0 ! Const for multiplying 
      if ( L > 0 ) then
          write(*,*) dr
          do b=1, nb
            r_lo = REAL( b-1 ) * dr
            r_hi = r_lo + dr
            h_id(b) = const * ( r_hi**3.0 - r_lo ** 3.0 ) ! Ideal Number
          enddo
      else
          do i=1,rconfigs
              do b=1, nb
                r_lo = real( b-1 ) *nptdr(i)
                r_hi = r_lo + nptdr(i)
                npth_id(i,b) = const * (r_hi**3.0 - r_lo ** 3.0)
              enddo
          enddo
      endif

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
            if (eweight == 1) eavg = eavg + en(cnt)
            do j=1, startskip
                read(12,*)
            enddo
            do j=1,natoms
                read(12,*) ctmp, (rtmp(a), a=1,3)
                ! Creates arrays of coordinates (in box coords)
                if (typ(j) == selec1) then
                    if ( L > 0.0 ) then 
                        r1(cnt,cnt1,:) = rtmp(:)/L
                    else
                        r1(cnt,cnt1,:) = rtmp(:)/nptL(cnt)
                    endif
                    mol1(cnt1) = mol(j)
                    ! Checks if selec1 and selec2 are the same
                    if(selec2 == selec1) then
                        if ( L > 0.0 ) then
                            r2(cnt,cnt1,:) = rtmp(:)/L
                        else
                            r2(cnt,cnt1,:) = rtmp(:)/nptL(cnt)
                        endif
                        mol2(cnt1) = mol(j)
                    endif
                    cnt1 = cnt1 + 1
                ! Note - this will never happen if selec1 and selec2 are
                ! equivalent, as that is handeled previously.
                else if (typ(j) == selec2) then
                    if ( L > 0.0 ) then
                        r2(cnt,cnt2,:) = rtmp(:)/L
                    else
                        r2(cnt,cnt2,:) = rtmp(:)/nptL(cnt)
                    endif
                    mol2(cnt2) = mol(j)
                    cnt2 = cnt2 + 1
                end if
            enddo
        else
            if (eweight == 1) read(13,*)
            do j=1, startskip
                read(12,*)
            enddo
            do j=1, natoms
                read(12,*)
            enddo
        end if
      enddo
      close(12)
      ! Calculate fluctuations in energy and blocking
      if (eweight == 1) then
          close(13)
          eavg = eavg/real(rconfigs)
          write(*,*) 'total average energy is ',eavg
          desqavg=0.0; den=0.0; decubavg=0.0
          do i=1,rconfigs
              den(i) = (en(i)-eavg)
              desqavg = desqavg + (den(i)**2)/real(rconfigs)
              decubavg= decubavg + (den(i)**3)/real(rconfigs)
          enddo
          write(*,*) 'total average fluctuation in sq energy is', desqavg
          write(*,*) 'total average fluctuation in cub energy is', decubavg
          eavgbl=0.0
          desqavgbl=0.0; decubavgbl=0.0
          do j=1,nblocks
            do i=1,rconfigs/nblocks
                k=(j-1)*rconfigs/nblocks+i
                eavgbl(j) = eavgbl(j) + en(k)
            end do
            eavgbl(j) = eavgbl(j)/(real(rconfigs)/real(nblocks))
            do i=1,rconfigs/nblocks
                k=(j-1)*rconfigs/nblocks+i
                desqavgbl(j) = desqavgbl(j) + (en(k)-eavgbl(j))**2
                decubavgbl(j)= decubavgbl(j) + (en(k)-eavgbl(j))**3
            end do
            desqavgbl(j) = desqavgbl(j)/(rconfigs/nblocks)
            decubavgbl(j) = decubavgbl(j)/(rconfigs/nblocks)
            write(*,*) 'block ',j,' average energy is ', eavgbl(j)
            write(*,*) 'block ',j,' average sq. fluct is ', desqavgbl(j)
            write(*,*) 'block ',j,' average cub. fluct is ', decubavgbl(j)
          enddo
      endif
      write(*,*) "Calculating GofR"
      ! Loop over configurations
      !$OMP PARALLEL DO schedule(static) DEFAULT(NONE) &
      !$OMP& PRIVATE(i, j, k, rij, rij_sq, b, nb, n1, n2, dr, cnt) &
      !$OMP& SHARED(nconfigs, selec1, selec2, r1, r2, g, h, h_id,rconfigs) 
      open(50,file="g_vals.dat")
      open(51,file="dr_vals.dat")
      do i=1, rconfigs
        do j=1, n1
            if (selec1 /= selec2) then
                do k=1, n2
                    if (mol1(j) /= mol2(k)) then
                        rij(:)      = r1(i,j,:) - r2(i,k,:)
                        rij(:)      = rij(:) - ANINT( rij(:) )
                        rij_sq      = SUM( rij**2 )
                        if ( L > 0 ) then
                            b           = FLOOR( SQRT( rij_sq ) / dr ) + 1
                        else
                            b           = FLOOR( SQRT( rij_sq ) / nptdr(i) ) + 1
                        endif
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
                            if ( L > 0 ) then
                                b       = FLOOR( SQRT( rij_sq ) / dr ) + 1
                            else
                                b       = FLOOR( SQRT( rij_sq ) / nptdr(i) ) + 1
                            endif
                            IF ( b <= nb ) h(i,b) = h(i,b) + 2
                        endif
                    enddo
                end if
            end if
        enddo
        ! Final Calculations
        g(i,:) = REAL( h(i,:) ) / REAL( n1 ) ! Norm by number of n1
        if ( L > 0.0 ) then
            g(i,:) = g(i,:) / h_id(:) ! Normalize based on ideal result
        else 
            g(i,:) = g(i,:) / npth_id(i,:) ! Normalize based on ideal result
        endif
        do b=1,nb
            write(50,*) g(i,b)
        enddo
        if (eweight == 1) eg(i,:) = en(i) * g(i,:) ! Weight by energies
        if (eweight == 1) e2g(i,:) = en(i)**2*g(i,:) ! Weight by sq energies
        if (eweight == 1) e3g(i,:) = en(i)**3*g(i,:) ! Weight by cub energies
      enddo
      do b=1,nb
        write(51,*) (REAL(b)-0.5)*dr*L
      enddo
      !$OMP END PARALLEL DO
      write(*,*) "Tidying things up!"
      write(censtr, "(I0)") selec1
      write(outstr, "(I0)") selec2
      ! Total calculation
      sumweight=0.0; weight=0.0 ! zeros reweighting
      do i=1, rconfigs
        if (reweight == 1) then
            weight = exp(-(en(i)-eavg)/(2*desqavg))
            sumweight = sumweight + weight
        endif
        do b=1,nb
            if (reweight == 0) then
                tot_gofr(b) = tot_gofr(b) + g(i,b)/rconfigs
                if (eweight == 1) tot_egofr(b)  = tot_egofr(b) + (eg(i,b)-eavg*g(i,b))/real(rconfigs)
                if (eweight == 1) tot_e2gofr(b) = tot_e2gofr(b) + ((den(i)**2.-desqavg)*g(i,b))/real(rconfigs)
                if (eweight == 1) tot_e3gofr(b) = tot_e3gofr(b) &
                & + ((den(i)**3. - decubavg -3*desqavg*den(i))*g(i,b))/real(rconfigs)
            else if (reweight == 1) then
                tot_gofr(b) = tot_gofr(b) + g(i,b)*weight
                if (eweight == 1) tot_egofr(b)  = tot_egofr(b) + (eg(i,b)-eavg*g(i,b))*weight
                if (eweight == 1) tot_e2gofr(b) = tot_e2gofr(b) + ((den(i)**2.-desqavg)*g(i,b))*weight 
                if (eweight == 1) tot_e3gofr(b) = tot_e3gofr(b) + ((den(i)**3. - decubavg -3*desqavg*den(i))*g(i,b))/real(rconfigs)*weight
            endif
        enddo
      enddo
      if (reweight == 1) then
          tot_gofr(:) = tot_gofr(:)/sumweight
          tot_egofr(:) = tot_egofr(:)/sumweight
          tot_e2gofr(:) = tot_e2gofr(:)/sumweight
          tot_e3gofr(:) = tot_e3gofr(:)/sumweight
      endif
      !tot_gofr(:) =dd tot_gofr(:)/rconfigs 
      !if (eweight == 1) tot_egofr(:) = tot_egofr(:)/rconfigs - eavg*tot_gofr(:)
      ! Block calculation
      write(*,*) "Jiggying up blocks!"
      if ( L > 0.0 ) dr = dr * L
      do j=1, nblocks
        bl_gofr=0.0
        if (eweight == 1) then
            bl_egofr=0.0; bl_e2gofr=0.0; bl_e3gofr=0.0
        endif
        do i=1,rconfigs/nblocks
            k=(j-1)*rconfigs/nblocks+i
            do b=1,nb
                bl_gofr(b) = bl_gofr(b) + g(k,b)/(rconfigs/nblocks)
                if (eweight == 1) bl_egofr(b) = bl_egofr(b) +(eg(k,b)-eavgbl(j)*g(k,b))/(rconfigs/nblocks)
                if (eweight == 1) bl_e2gofr(b) = bl_e2gofr(b) +((den(k)**2.-desqavgbl(j))*g(k,b))/(rconfigs/nblocks)
                if (eweight == 1) bl_e3gofr(b) = bl_e3gofr(b) + ((den(k)**3. -decubavgbl(j) -3*desqavgbl(j)*den(k))*g(k,b))/real(rconfigs/nblocks)
            enddo
        enddo
        
        write(blstr, "(I0)") j
        open(20+j,file='bl_'//trim(blstr)//'_pairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
        if (eweight == 1) then
            open(42+j,file='bl_'//trim(blstr)//'_'//trim(etype)//'pairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
        endif
        do b=1,nb
            write(20+j,'(2f15.8)' ) (REAL(b)-0.5)*dr, bl_gofr(b)
            if (eweight == 1) then
                write(42+j,'(4f15.8)' ) (REAL(b)-0.5)*dr, bl_egofr(b), bl_e2gofr(b), bl_e3gofr(b)
            endif
        enddo
        close(20+j)
        if (eweight == 1) close(42+j)
      enddo


      
      write(*,*) "Converted back to original coordinates"
      ! Output to file
      write( *,* ) 'Outputing to output file'
      open(20, file='pairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
      if (eweight == 1) open(42,file=trim(etype)//'pairdist_'//trim(censtr)//'_'//trim(outstr)//'.dat')
      do b=1, nb
        write(20,'(1f15.8,5x,1e15.6)' ) (REAL(b)-0.5)*dr, tot_gofr(b)
        if (eweight == 1) write(42,'(1f15.8,5x,1e15.6,1e15.6, 1e15.6)' ) (REAL(b)-0.5)*dr, tot_egofr(b), tot_e2gofr(b), tot_e3gofr(b)
      enddo
      close(20)
      if (eweight == 1) close(42)
      write(*,*) "GofR Calculation Complete"
      if (eweight == 1) write(*,*) "Weighted GofR Calculation Complete"
     


END Program pairdist


