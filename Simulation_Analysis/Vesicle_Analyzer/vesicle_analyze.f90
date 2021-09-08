module constants
  implicit none

  real, parameter :: pi=DACOS(-1.D0)


end module constants

module parameters
  integer :: lipstart, lipstop,frame
  integer :: ufile, lfile, drhfile,drfile, tfile, hdrfile, ldrfile
  integer :: nmol, nframes, atoms_per_mol, hg_index, lg_index, natoms
  real :: L
  real :: rmin, rmax, rbin
  character(len=20) :: lengthfile, trajfile
  common lipstart, lipstop, L, natoms, hg_index, lg_index
  common atoms_per_mol, nmol, lengthfile, nframes
  common tfile, lfile, drfile, trajfile
  common rmin, rmax, rbin
end module parameters

Program CalcVesc
  Use, INTRINSIC :: iso_fortran_env, ONLY: input_unit, output_unit, &
    error_unit, iostat_end, iostat_eor
  use constants
  use parameters
  implicit none
  integer :: ftype
  integer :: i, j, k
  real, dimension(3) :: comr
  real, allocatable :: r(:,:)
  call Read_Input()
  ! Allocations
  allocate(r(natoms,3))
  ! End Allocations

  !Input
  tfile=11; lfile=12; 
  open(tfile, file=trim(trajfile), status='old')
  open(lfile, file=trim(lengthfile), status='old')
  ! Output
  drfile=21; drhfile=22
  hdrfile=23; ldrfile=24
  open(drhfile,file="drh_out.dat")
  open(drfile, file="dr_out_hist_raw.dat")
  open(hdrfile, file="hdr_out_hist_raw.dat")
  open(ldrfile, file="ldr_out_hist_raw.dat")
  ! Operations on that frame
  do frame=1,nframes
    call ReadFrame(11,12,2,r)
    call Calc_VescAltCOM(r,comr)
    call Calc_RadialProfile(r,comr)  
  enddo
  close(tfile)
  close(lfile)
  close(drhfile)
  close(drfile)
  write(*,*) "End Program"
End Program CalcVesc

Subroutine Read_Input()
  use parameters
  implicit none
  character (len=6) :: ctmp
  character(len=20) :: sysfile
  open(10,file="vesicle.in", status='old')
  ! Read Input Parameters
  read(10,*)
  read(10,*) sysfile, trajfile,lengthfile
  read(10,*)
  read(10,*) nmol, atoms_per_mol
  read(10,*)
  read(10,*) hg_index, lg_index
  read(10,*)
  read(10,*) nframes
  read(10,*)
  read(10,*) rmin, rmax, rbin
  lipstart=1
  lipstop=nmol*atoms_per_mol
  close(10)
  open(10, file=trim(sysfile), status='old')
  read(10,*)
  read(10,*) natoms, ctmp
  close(10)
End Subroutine Read_Input

Subroutine Calc_RadialProfile(r,comr)
  use constants
  use parameters
  implicit none
  integer :: i, j, atom,a_id, h_id,l_id
  integer :: cnt_inner, cnt_outer
  real :: M, dr_outer, dr_inner
  real, dimension(3) :: comr, drtmp_h, drtmp_l,drtmp
  real, dimension(nmol) :: dr_h, dr_l
  real, dimension(natoms) :: drsq, dr 
  real, dimension(natoms,3) :: r
  integer, dimension(nmol*atoms_per_mol) :: rprofile
  drtmp_h = 0.0; drtmp_l=0.0
  dr_h = 0.0; dr_l = 0.0
  cnt_inner = 0; cnt_outer = 0
  dr_outer = 0.0; dr_inner = 0.0
  drsq = 0.0; dr =0.0
  h_id=1; l_id = 1
  do i=1,nmol
    h_id = (i-1)*atoms_per_mol+hg_index
    l_id = (i-1)*atoms_per_mol+lg_index
    do atom=1,atoms_per_mol
      a_id = (i-1)*atoms_per_mol + atom
      do j=1,3
        drtmp(j) = r(a_id,j) -comr(j)
        drtmp(j) = drtmp(j) + L*anint(drtmp(j)/L)
        drsq(a_id) = drsq(a_id) + drtmp(j)**2.
        dr(a_id) = sqrt(drsq(a_id))
      enddo !j
      if (h_id .eq. a_id) then
        dr_h(i) = dr(a_id)
      endif
      if (l_id .eq. a_id) then
        dr_l(i) = dr(a_id)
      endif
    enddo !atom
    if (dr_h(i) > dr_l(i)) then
      dr_outer = dr_outer + sqrt(dr_h(i))
      cnt_outer = cnt_outer + 1
    else
      dr_inner = dr_inner + sqrt(dr_h(i))
      cnt_inner = cnt_inner + 1
    endif
  enddo!i
  call Histogram(dr, nmol*atoms_per_mol, rmin, rmax, rbin, drfile)
  call Histogram(dr_h, nmol, rmin, rmax, rbin, hdrfile)
  call Histogram(dr_l, nmol, rmin, rmax, rbin, ldrfile)
  write(drhfile,*) frame, dr_outer/real(cnt_outer), dr_inner/real(cnt_inner)
End Subroutine Calc_RadialProfile

Subroutine Histogram(values,n,min_val, max_val, nbins,fout)
  integer :: i, j, n, cnt, bin, fout
  real  :: min_val, max_val, nbins,bwidth
  real, dimension(n) :: values
  integer, dimension(nbins) :: hist
  bwidth = (max_val-min_val)/nbins
  hist=0; cnt=0
  do i=1, n
    bin = ceiling((values(i)-min_val)/bwidth)
    if (bin > nbins) then
      stop "Error: Bins should be widened"
    endif
    if (bin .eq. 1) then
      write(*,*) "first", values(i)
    endif
    hist(bin) = hist(bin) + 1
    cnt = cnt + 1
  enddo !i
  do i=1, nbins
    write(fout,*) i*bwidth-bwidth/2.0+min_val, hist(i)
  enddo!i

End Subroutine

Subroutine Calc_VescCOM(r, comr)
  use constants
  use parameters
  implicit none
  integer :: i, j
  real :: M
  real, dimension(3) :: comr
  real, dimension(natoms,3) :: r
  comr = 0.0
  M = 0.0
  do i=lipstart,lipstop
    do j = 1, 3
      comr(j) = comr(j) + r(i,j)*72.0
    enddo !j
    M = M + 72.0
  enddo !i

End Subroutine

Subroutine Calc_VescAltCOM(r,comr)
  ! This is a subroutine for calculating the center of mass
  ! regardless of periodic boundary conditions, as is described here
  ! https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
  ! Note that this works by transforming each dim into polar coords
  ! and then transforming back. 
  ! Key output is comr, which stores the center of mass location.
  use constants
  use parameters
  implicit none
  integer :: i, j
  real :: M
  real, dimension(3) :: comr,thta,a1,a2,thtav
  real, dimension(natoms,3) :: r
  comr = 0.0; thta = 0.0
  M = 0.0
  a1 = 0.0; a2 = 0.0; thtav=0.0
  do i=lipstart, lipstop
    do j=1,3
      thta(j) = r(i,j)/L*2*pi
      a1(j) = a1(j) +  72.0*cos(thta(j))
      a2(j) = a2(j) +  72.0*sin(thta(j))
    enddo !j
    M = M + 72.0
  enddo !i
  do j=1,3
    a1(j) = a1(j)/M
    a2(j) = a2(j)/M
    thtav(j) = atan2(-a2(j),-a1(j)) + pi
    comr(j) = L*thtav(j)/(2*pi)
    if (comr(j) > L/2) then
      comr(j) = comr(j)-L
    endif
  enddo !j
End Subroutine

Subroutine  ReadFrame(ufile,lfile,ftype,r)
  ! This reads xyz or lammps dump files (of type id type xu yu zu)
  ! outputs an array of size (natoms, 3) with all coordinates
  use constants
  use parameters
  implicit none
  integer :: i, j, k
  integer :: ufile, ftype, lfile
  real :: ctmp
  real, dimension(natoms,3) :: r
  r=0.0
  if (mod(frame,100) .eq. 0) then
    write(*,*) frame
  endif
  read(lfile,*) L
  if (ftype .eq. 1) then
    do i = 1, 2
      read(ufile,*)
    enddo
    do i = 1, natoms
        read(ufile,*) ctmp, r(i,1), r(i,2), r(i,3)
    enddo ! i
  else if (ftype .eq. 2) then
    do i = 1, 9
      read(ufile,*)
    enddo
    do i = 1, natoms
    read(ufile,*) ctmp, ctmp, r(i,1), r(i,2), r(i,3)
    enddo ! i
  else
    write(*,*) "Error: ftype provided is not [1] xyz or [2] lammps dump"
    stop
  endif
End Subroutine ReadFrame
