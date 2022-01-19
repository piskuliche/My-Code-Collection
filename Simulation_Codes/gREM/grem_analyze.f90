! Copyright 2022 Ezekiel Piskulich, Boston University
!
! This software is used to analyze a single walker from a grem simulation after having
! sorted using the grem_sort.f90 code.
! 
! This code takes a single input file (../grem.input) in the directory above from where the code has been run.
! an example of this input file is included in the examples directory.
!
! For this code to work - walker files should be labeled as walker-XXX.lammpsdump
! Output files: P2-XXX.dat area-XXX.dat thick-XXX.dat
!
! For COM calculations, assumes MARTINI standard of 72 g/mol 


module constants
    implicit none
    real, parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: numfile
  integer :: fstart, fstop, nwindows,natoms, nlipids
  integer :: nlog, dumpfreq, ntailpair
  integer :: hatom, lt1, lt2, lt3,lt4
  real :: xlo, ylo, zlo, xhi, yhi, zhi
  real, dimension(3) :: L
  real, dimension(50) :: ltail
  common fstart, fstop, nwindows,natoms,nlog,dumpfreq,nlipids
  common L, hatom, ntailpair, ltail
  common xlo, ylo, zlo, xhi, yhi, zhi
end module parameters

character(len=20) function str(z)
  !   "Convert an integer to string."
  integer, intent(in) :: z
  write (str, *) z
  str = adjustl(str)
end function str

program analyze_grem
  use constants
  use parameters
  implicit none
  integer :: i, j, k
  integer :: outvalue, nfile
  integer :: nframes
  integer, dimension(200) :: sortdata, infiles, outfiles

  character(len=20) :: str,ctmp

  infiles=0
  outfiles=0
  read(*,*) nfile
  write(*,*) "Begin Program"
  call read_input()
  write(*,*) fstart,fstop,nwindows,nlog

  open(100, file="walker-"//trim(str(nfile))//".lammpsdump",status='old')
  open(101, file="area-"//trim(str(nfile))//".dat")
  open(102, file="thick-"//trim(str(nfile))//".dat")
  open(103, file="p2-"//trim(str(nfile))//".dat")
 
  nframes=nlog/dumpfreq
  do i=fstart,fstop
    do j=1,nframes
      call readframe(100)
    enddo !j
  enddo !i

  close(100)
  close(101)
  close(102)
  close(103)

end program analyze_grem

subroutine readframe(infile)
  use constants
  use parameters
  implicit none
  integer :: i,j,k
  integer :: atom1, atom2
  integer :: infile
  integer :: areafile, p2file, thickfile
  integer :: id, mol, atyp
  integer,dimension(3) :: ival

  real :: area, dz, comza, comzb
  real :: thick, drsq1, drsq2
  real :: Ma, Mb, M
  real :: p2, p2tmp
  real, dimension(3) :: dr, drtmp,defz
  real, dimension(50,3) :: r

  character(len=200) :: ctmp

  areafile = infile+1
  thickfile= infile+2
  p2file = infile+3
  r=0.0; ival=0
  Ma=0.0; Mb=0.0
  thick=0.0; comza=0.0; comzb=0.0
  defz = 0; p2 = 0.0
  defz(3) = 1
  ! header
  do i=1,9
    read(infile,'(A)') ctmp
    if (i==6) read(ctmp,*) xlo, xhi
    if (i==7) read(ctmp,*) ylo, yhi
    if (i==8) read(ctmp,*) zlo, zhi
  enddo !i
  L(1) = xhi-xlo; L(2) = yhi-ylo; L(3)=zhi-zlo
  area = L(1)*L(2)
  write(areafile,*) area
  do i=1,nlipids
    do j=1,natoms
      read(infile,'(A)') ctmp
      read(ctmp,*) id, mol, atyp, r(j,1), r(j,2), r(j,3), ival(1), ival(2), ival(3)
      do k=1,3
        r(j,k) = r(j,k) + ival(k)*L(k)
      enddo
    enddo !j
    call pbc_dist(r(hatom,:), r(lt1,:),L,dr)
    ! Check whether pointing up or down
    if (dr(3) > 0) then
      Ma = Ma + 72.0
      comza = comza +  r(hatom,3) * 72.0
    else
      Mb = Mb + 72.0
      comzb =comzb + r(hatom,3) * 72.0
    endif
    do k=1,ntailpair
      atom1 = ltail(2*(k-1)+1)
      atom2 = ltail(2*k)
      call pbc_dist(r(atom1,:),r(atom2,:),L,drtmp)
      call calc_p2(drtmp,defz,p2tmp)
      p2 = p2 + p2tmp
    enddo
  enddo !i
  ! Calculate Thickness
  comza = comza/Ma
  comzb = comzb/Mb
  dz = comza - comzb
  dz = dz - L(3)*anint(dz/L(3))
  thick = thick + dz
  ! Calculate P2
  p2 = p2/real(nlipids*ntailpair)
  write(thickfile,*) abs(thick)
  write(p2file,*) p2
end subroutine readframe

subroutine calc_p2(dr,defz,p2)
  implicit none
  integer :: i
  real :: drsq
  real :: costhta
  real :: p2
  real,dimension(3) :: dr, defz

  drsq = dot_product(dr,dr)
  do i=1,3
    dr(i) = dr(i)/sqrt(drsq)
  enddo
  costhta = dot_product(dr,defz)
  p2 = (3*costhta**2. -1)/2.0
end subroutine calc_p2


subroutine pbc_dist(r1,r2,L,dr)
  implicit none
  integer :: i
  real, dimension(3) :: r1, r2, L, dr
  dr = 0.0
  do i=1,3
    dr(i) = r1(i) - r2(i)
    dr(i) = dr(i) - L(i)*anint(dr(i)/L(i))
  enddo
end subroutine pbc_dist


subroutine readwalkers(sortdata)
  use constants
  use parameters
  implicit none
  integer :: t,i,j,k
  integer, dimension(200) :: sortdata
  sortdata=0
  read(12,*) t, (sortdata(i), i=1,nwindows)
end subroutine

subroutine read_input()
  use parameters
  implicit none
  
  integer :: infile,i
  infile = 11
  ltail=0
  open(infile, file="../grem.input", status='old')
  read(infile,*)
  read(infile,*) nwindows
  read(infile,*)
  read(infile,*) fstart, fstop
  read(infile,*)
  read(infile,*) nlog, dumpfreq
  read(infile,*) 
  read(infile,*) natoms, nlipids
  read(infile,*) 
  read(infile,*) hatom, ntailpair
  read(infile,*)
  read(infile,*) (ltail(i), i=1,ntailpair*2)
  close(infile)
end subroutine read_input
