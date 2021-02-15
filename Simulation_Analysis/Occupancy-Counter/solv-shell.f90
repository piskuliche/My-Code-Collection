Program Solv_Finder
      implicit none
      integer :: i, j, k, ntimes, t, sk,cnt
      integer :: acry_natms, nwater, nacryl
      real :: volume,maxdrsq
      integer, dimension(100) :: h2o_count, ad,h1, h2, hbond_count
      real, dimension (100) :: criteria,critsq
      character(len=3), dimension(100) ::  at_name
      character(len=2) :: ctmp2
      character(len=3) :: ctmp, ctmp3
      character(len=20) :: nfile
      real, dimension(3) :: L
      real, dimension(1000,3) :: rO, r1, r2
      real, dimension(100,100,3) :: racryl
      
      
      ! Set atomkey to 0.0  
      criteria=0.0
      rO=0.0; r1 = 0.0; r2 = 0.0
      racryl=0.0
      h2o_count=0
      critsq=0.0
      criteria=0.0
      h1 = 0; h2 = 0; ad = 0
        
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
        read(11,*) at_name(i), criteria(i), ad(i), h1(i), h2(i)
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
      open(14,file="hbonds.dat")
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
            read(12,*) ctmp, (r1(i,k),k=1,3)
            read(12,*) ctmp, (r2(i,k),k=1,3)
            r1(i,:) = r1(i,:) - L(:)*anint(( r1(i,:) - rO(i,:) )/L(:))
            r2(i,:) = r2(i,:) - L(:)*anint(( r2(i,:) - rO(i,:) )/L(:))
        enddo
        ! Finds the distance partners at that timestep
        call find_partners(nacryl, nwater, acry_natms,L, rO,r1,r2, racryl, critsq, h2o_count,h1,h2,ad,hbond_count)
        write(13,'(100I8.4)') (h2o_count(k),k=1,acry_natms)
        write(14,'(100I8.4)') (hbond_count(k), k=1,acry_natms)
      enddo
      close(12)
      close(13)
      close(14)


End Program   

Subroutine find_partners(nacryl, nwater, acry_natms,L, rO,r1,r2, racryl, critsq, h2o_count,h1, h2,ad,hbond_count)
    implicit none
    integer :: i, j, k, cnt
    integer :: nacryl, nwater, acry_natms,hbonds
    integer, dimension(100) :: h1, h2,ad
    real :: drsq, maxdrsq
    real, dimension(3) :: L
    real, dimension(100,100,3) :: racryl
    real, dimension(1000,3) :: rO,r1,r2
    real, dimension(3) :: dr
    real, dimension(100) :: critsq
    real, dimension(500) :: h2o_mindist
    integer, dimension(500) :: h2o_loc
    integer, dimension(100) :: h2o_count, hbond_count
    integer, dimension(1000) :: counted
    maxdrsq = 3*(L(1)/2)**2.0

    h2o_count=0; hbond_count=0
    counted=0
    h2o_loc = 0
    hbonds=0
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
                    ! Checks H-Bond Partners
                    if ( ad(i) .eq. 1 ) then
                        call hbond_partners(i,j,k,dr,L,rO,r1,r2,racryl,h1,h2,hbonds)
                        hbond_count(i) = hbond_count(i) + hbonds
                    endif
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

Subroutine hbond_partners(i,j,k,drOX,L,rO,r1,r2,racryl,h1,h2,hbonds)

    implicit none
    integer :: i, j, k, cnt
    integer :: hbonds
    integer :: nacryl, nwater, acry_natms
    integer, dimension(100) :: h1, h2
    real :: dOX, dHX1, dHX2, dHO1, dHO2, ang
    real :: rOXmax, rHXmax,angmax
    real, dimension(3) :: drOX, drHX1, drHX2,drHO1, drHO2, eox, ehx1, ehx2, eho1, eho2
    real, dimension(3) :: L
    real, dimension(100,100,3) :: racryl
    real, dimension(1000,3) :: rO, r1, r2
    ! Zero stuff
    drHX1 = 0.0; drHX2 = 0.0; eox = 0.0; ehx1 = 0.0; ehx2 = 0.0
    ! H-Bond Criteria
    rOXmax = 3.5; rHXmax = 2.45; angmax=30
    ! Note i = acryl atom, j = acryl mol, k = water molecule
    ! Zeros hbonds count
    hbonds=0
    ! Calculates the OX distance
    dOX = sqrt(dot_product(drOX,drOX))
    if  ( dOX < rOXmax ) then
        ! Checks Water Donations to Acryl
        drHX1(:)=r1(k,:)-racryl(i,j,:) - L(:)*anint((r1(k,:)-racryl(i,j,:))/L(:))
        drHX2(:)=r2(k,:)-racryl(i,j,:) - L(:)*anint((r2(k,:)-racryl(i,j,:))/L(:))
        ! Water H, Acryl X (O or N) distance
        dHX1=sqrt(dot_product(drHX1,drHX1))
        dHX2=sqrt(dot_product(drHX2,drHX2))
        eox(:) = drOX(:)/dOX
        if ( dHX1 .lt. rHXmax ) then
            ehx1(:) = drHX1(:)/dHX1
            ang = acosd(dot_product(eox,ehx1))
            if ( ang < angmax ) then
                ! H-bonded
                hbonds = hbonds + 1
            endif
        endif
        if ( dHX2 .lt. rHXmax ) then
            ehx2(:) = drHX2(:)/dHX2
            ang = acosd(dot_product(eox,ehx2))
            if ( ang < angmax ) then
                ! H-bonded
                hbonds = hbonds + 1 
            endif
        endif
        ! Checks Acryl Donations to Water
        if ( h1(i) .ne. 0) then
            drHO1(:)=racryl(h1(i),j,:) - rO(k,:) - L(:)*anint((racryl(h1(i),j,:) -rO(k,:))/L(:))
            ! Water H, Acryl X (O or N) distance
            dHO1=sqrt(dot_product(drHO1,drHO1))
            if ( dHO1 .lt. rHXmax ) then
                eho1(:) = drHO1(:)/dHO1
                ang = acosd(dot_product(-eox,eho1))
                if ( ang < angmax ) then
                    hbonds = hbonds + 1
                endif
            endif
        endif
        if ( h2(i) .ne. 0) then
            drHO2(:)=racryl(h2(i),j,:) - rO(k,:) - L(:)*anint((racryl(h2(i),j,:)-rO(k,:))/L(:))
            dHO2=sqrt(dot_product(drHO2,drHO2))
            if ( dHO2 .lt. rHXmax ) then
                eho2(:) = drHO2(:)/dHO2
                ang = acosd(dot_product(-eox,eho2))
                if ( ang < angmax ) then
                    hbonds = hbonds + 1
                endif
            endif
        endif
    endif

    

End Subroutine
