! This code does a radius of gyration calculation for a collection of atoms

Program calc_gyro
   Use, INTRINSIC :: iso_fortran_env, ONLY: input_unit, output_unit, &
           error_unit, iostat_end, iostat_eor
   Implicit None 
   !**********
   !*Integers*
   !**********
   ! Loop Variables
   integer :: i, j, k, a 
   ! Counters
   integer :: n, cnt, cntn
   ! Input Values
   integer :: startconfig, endconfig, startskip, endskip
   integer :: selec, nblocks, natoms, bstart, blen
   integer :: nconfigs, or_int, rconfigs
   integer :: ioerr

   real :: ctmp, totRg, drsq

   ! Arrays
   real, dimension(3) :: rij, rcm, rtmp
   real, dimension(3) :: L
   real, allocatable :: r(:,:,:)
   real, allocatable :: Rg(:), typ(:), mol(:), Rgbl(:)
   


   ! Filenames
   character(len=20) nfile, molfile

   ! Used Files
   ! 11: Open Trajectory File

   NAMELIST /nml/ nfile, molfile, startconfig, endconfig, startskip, endskip,&
       & selec, natoms, nblocks

   ! Example Defaults
   nfile = "traj.xyz"
   molfile = "molfile.info"
   startconfig = 0
   endconfig   = 10
   startskip   = 2
   endskip     = 0
   selec       = 2 
   natoms      = 100


   ! Reads the namelist from standard input
   read ( unit=input_unit, nml=nml, iostat=ioerr)
   if ( ioerr /= 0 ) then
     write( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from stdin', ioerr
     if ( ioerr == iostat_eor ) write( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write( unit=error_unit, fmt='(a)') 'End of file'
   end if

   write(*,*) '*******START*******'
   write(*,*) 'nfile: ', nfile
   write(*,*) 'molfile: ', molfile
   write(*,*) 'startconfig: ', startconfig
   write(*,*) 'endconfig: ', endconfig
   write(*,*) 'startskip: ', startskip
   write(*,*) 'endskip: ', endskip
   write(*,*) 'selec: ', selec
   write(*,*) 'natoms: ', natoms
   write(*,*) 'nblocks: ', nblocks
   
   ! This section reads in the molecule and the type from the
   ! molfile created by the setup script.
   ! It also counts the total number that match the selection (n)
   allocate(typ(natoms), mol(natoms))
   rconfigs = startconfig-endconfig
   allocate(Rg(rconfigs))
   allocate(Rgbl(nblocks))

   mol=0; typ=0; n = 0
   open(11, file=trim(molfile), status='old')
   do j=1, natoms
       read(11,*) ctmp, mol(j), typ(j), ctmp
       if (typ(j) == selec) then
           n = n + 1
       endif
   enddo
   close(11)
   write(*,*) "There are ", n, " of your selection"
   
   ! Last allocation, r
   allocate(r(rconfigs,n,3)) 

   open(12, file=nfile, status="old")
   do i=1, startconfig
       do j=1, startskip
           read(12,*)
       end do
       do j=1, natoms
           read(12,*)
       end do
   end do
   cnt = 0
   do i=1, nconfigs
    if (MOD(i,1000) == 0) then
        write(*,*) 'Configuration ',i,' read in.'
    end if
    if (MOD(i,or_int) == 0) then
        cnt = cnt + 1
        cntn = 1
        do j=1,startskip
            read(12,*)
        enddo
        do j=1,natoms
            read(12,*) ctmp, (rtmp(a),a=1,3)
            if (typ(j) == selec) then
                r(cnt, cntn, :) = rtmp(:)/L ! Convert to box coords
                if (cntn > 1) then ! Wrapping around the first atom in the frame
                    rij(:) = r(cnt,1,:) - r(cnt,cntn,:)
                    r(cnt,cntn,:) = r(cnt,cntn,:) - ANINT( rij(:) ) 
                endif
                cntn = cntn + 1             ! Increment counter
            endif
        enddo
    end if
   enddo
   close(12)
   write(*,*) "Calculating atomic Radius of Gyration"
   write(*,*) "There are ", n, " ion pairs in the system"
   Rg = 0
   ! Loops over configurations
   do i=1, rconfigs
       rij=0; rcm=0
       ! Calculates Center of Mass of Ion Cluster
       ! Note I still need this, otherwise I won't be looking from center of 
       ! the ion cluster
       do j=1, n
           rcm(:) = rcm(:) + r(i,j,:)
       enddo
       rcm(:) = rcm(:) / real(cntn)
       ! Loops over ions
       do j=1, n
           rij(:) = r(i,j,:) - rcm(:) ! Subtracts the center of mass position
           drsq   = rij(1)**2 + rij(2)**2 + rij(3)**2
           Rg(i) = Rg(i) + drsq
       enddo
       Rg(i) = Rg(i)/real(n)
   enddo
   ! Calculate Total and block average
   totRg = SUM( Rg ) / rconfigs
   blen = rconfigs/nblocks
   do i=1,nblocks
       bstart = (i-1)*blen
       do j=1, blen
        Rgbl(i) = Rgbl(i) + Rg(bstart+j)
       enddo
       Rgbl(i) = Rgbl(i)/blen
   enddo
   
   ! Write out final results
   write(*,*) totRg, (Rgbl(k), k=1,nblocks)


end program calc_gyro

    
    
    
