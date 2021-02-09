Program Solv_Finder
      implicit none
      integer :: i, j, k, ntimes, t, sk,cnt
      integer :: acry_natms, nwater, nacryl
      real :: volume,maxdrsq
      integer, dimension(100) :: h2o_count
      real, dimension (100) :: criteria,critsq
      character(len=3), dimension(100) ::  at_name
      character(len=2) :: ctmp2
      character(len=3) :: ctmp, ctmp3
      character(len=20) :: nfile
      real, dimension(3) :: L
      real, dimension(1000,3) :: rO
      real, dimension(100,100,3) :: racryl
      
      
      ! Set atomkey to 0.0  
      criteria=0.0
      rO=0.0
      racryl=0.0
      h2o_count=0
      critsq=0.0
      criteria=0.0
        
      ! Reads the input file
      open(10,file='solv.in',status='old')
      read(10,*)
      read(10,*) nfile
      read(10,*) 
      read(10,*) volume, ntimes
      read(10,*)
      read(10,*) nwater, nacryl

      L(1)=volume ** (1.0/3.0)
      L(2)=L(1)
      L(3)=L(1)
      write(*,*) "Box length ", L(1)
      maxdrsq = 3*(L(1)/2)**2.0
      write(*,*) "Max Sq. Distance ", maxdrsq

      ! Reads the distances
      open(11,file="inp_vals", status='old')
      read(11,*) acry_natms
      write(*,*) "There are ", acry_natms, " per Acrylamide"
      cnt = 0
      do i=1,acry_natms
        read(11,*) at_name(i), criteria(i)
        critsq(i)=criteria(i)**2.0
        if (criteria(i) .ne. 0.0) then
            cnt = cnt + 1
        endif
      enddo
      close(11)
      write(*,*) "There are ", cnt," atom criteria to compare"

      ! Opens the trajectory
      open(12,file=nfile, status='old')
      open(13,file="solv_data.dat")
      do t=1,ntimes
        ! Skip Header
        do sk=1,9
            read(12,*)
        enddo
        do i=1,nacryl
            do j=1,acry_natms
                read(12,*) ctmp, (racryl(j,i,k),k=1,3)
            enddo
        enddo
        do i=1,nwater
            ! Read Oxygen Atoms
            read(12,*) ctmp, (rO(i,k),k=1,3)
            ! Skip H Atoms (for now)
            read(12,*)
            read(12,*)
        enddo
        ! Finds the distance partners at that timestep
        call find_partners(nacryl, nwater, acry_natms,L, rO, racryl, critsq, h2o_count)
        write(13,'(100I8.4)') (h2o_count(k),k=1,acry_natms)
      enddo
      close(12)
      close(13)


End Program   

Subroutine find_partners(nacryl, nwater, acry_natms,L, rO, racryl, critsq, h2o_count)
    implicit none
    integer :: i, j, k, cnt
    integer :: nacryl, nwater, acry_natms
    real :: drsq, maxdrsq
    real, dimension(3) :: L
    real, dimension(100,100,3) :: racryl
    real, dimension(1000,3) :: rO
    real, dimension(3) :: dr
    real, dimension(100) :: critsq
    real, dimension(500) :: h2o_mindist
    integer, dimension(500) :: h2o_loc
    integer, dimension(100) :: h2o_count
    integer, dimension(1000) :: counted
    maxdrsq = 3*(L(1)/2)**2.0

    h2o_count=0
    counted=0
    h2o_loc = 0
    h2o_mindist = 50.0
    ! loop over acry_natms
    do i=1,acry_natms
        ! loop over acry molecs
        if ( critsq(i) .ne. 0.0 ) then
            do j=1,nacryl
                ! loop over water molecs
                do k=1,nwater
                    ! calc distance
                    dr = rO(k,:)-racryl(i,j,:) - L(:)*anint((rO(k,:)-racryl(i,j,:))/L(:))
                    drsq = dr(1)**2.0 + dr(2)**2.0 + dr(3)**2.0
                    if (drsq .lt. critsq(i)) then
                        ! This double checks that it is he shortest distance the
                        ! h2o is involved in
                        if (drsq .lt. h2o_mindist(k)) then
                            h2o_count(i) = h2o_count(i) +  1
                            if (h2o_loc(k) .ne. 0) then
                                h2o_count(h2o_loc(k)) = h2o_count(h2o_loc(k))-1
                            endif
                            h2o_mindist(k) = drsq
                            h2o_loc(k) = i
                            counted(k)=1
                        endif
                    endif
                    if (drsq .gt. maxdrsq) then
                        write(*,*) drsq, maxdrsq, "Error"
                    endif
                enddo
            enddo
        endif
    enddo

    
End Subroutine
