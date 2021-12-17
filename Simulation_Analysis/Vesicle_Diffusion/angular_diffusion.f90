module constants
  implicit none
  real, parameter :: pi=DACOS(-1.D0)

end module constants

module parameters
  integer :: tfile
  integer :: natoms,nmol, atoms_per_mol, atom_index
  integer :: nframes, nstart, nstop, nsep, nskip
  real    :: L
  character(len=20) :: lengthfile, trajfile

  common tfile
  common nmol, atoms_per_mol, atom_index
  common nframes, nstart, nstop, nsep, nskip
  common L
  common lengthfile, trajfile

end module parameters


Program Angular
  use constants
  use parameters
  implicit none

  call Read_Input()
  write(*,*) "There are ", nframes, " frames"
End Program Angular



Subroutine ReadFrame()

  ! This reads xyz or lammps dump files (of type id type xu yu zu)
  ! outputs an array of size (natoms, 3) with all coordinates
  use constants
  use parameters
  implicit none
  integer :: i, j, k

  real :: dtmp, xlo, xhi
  real, dimension(natoms,3) :: r
  
  character(len=20) :: ctmp

  r=0.0
  do i=1,9
    read(tfile,*) ctmp
    if (i == 6) read(ctmp,*) xlo, xhi
  enddo ! i
  L = xhi - xlo
  do i = 1, natoms
    read(tfile,*) dtmp, dtmp, r(i,1), r(i,2), r(i,3)
  enddo ! i

End Subroutine ReadFrame

subroutine Read_Input()
  use parameters
  use constants
  implicit none
  character(len=6) :: ctmp
  character(len=20) :: sysfile
  open(10,file='angular.in', status='old')
  read(10,*)
  read(10,*) trajfile, lengthfile
  read(10,*)
  read(10,*) nmol, atoms_per_mol
  read(10,*) 
  read(10,*) atom_index
  read(10,*)
  read(10,*) nframes, nstart, nstop
  read(10,*)
  read(10,*) nsep, nskip
  close(10)
  open(11, file=trajfile, status='old')
  read(11,*)
  read(11,*) natoms
  close(11)

end subroutine Read_Input
