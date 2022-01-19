module constants
  implicit none
  real,parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: nspecies,nframes
  integer, dimension(10) :: atoms_per_species
  real, dimension(3) :: L
  integer :: atom_start, atom_stop, natoms
  integer :: ufile, lfile
  real :: ctmp, xlo, xhi, ylo, yhi, zlo, zhi

  common nspecies,nframes
  common atoms_per_species
  common L
  common atom_start, atom_stop, natoms
  common ufile
  common ctmp, xlo, xhi, ylo, yhi, zlo, zhi
end module parameters

Program Main
  use, intrinsic :: iso_fortran_env, ONLY: input_unit, output_unit, &
    error_unit, iostat_end, iostat_eor
  use constants
  use parameters
  implicit none
  integer i, j, k
  real, allocatable :: r(:,:)
  real, dimension(10,3) :: com_r

  atoms_per_species(:)=0 
  ufile=11; lfile = 12
  call Read_Input(atoms_per_species)
  allocate(r(natoms,3))
  r(:,:) = 0.0
  
  do i=1,nspecies
    write(*,*) atoms_per_species(i)
  end do

  do frame=1,nframes
    call Read_Frame(r)
    atom_start = 0; atom_stop = 0.0
    do i=1, nspecies
      atom_start = atom_stop + 1
      atom_stop = atom_start + atoms_per_species(i)
      call Calc_COM(r, atom_start, atom_stop, com_r(i,:))
    end do
  enddo
End Program

Subroutine Calc_COM(r,atom_start, atom_stop,tmpcom)
  ! This is a subroutine for calculating the center of mass regardless of PBC

  use constants 
  use parameters
  implicit none
  integer :: i, j
  integer :: atom_start, atom_stop
  real :: M
  real, dimension(natoms,3) :: r
  real, dimension(3) :: tmpcom, thta, a1, a2, thta_v
  tmpcom = 0.0; thta = 0.0
  M = 0.0 
  a1 = 0.0; a2 = 0.0; thta_v = 0.0
  do i=atom_start, atom_stop
    do j=1,3
      thta(j) = (r(i,j)*2*pi)/L(j)
      a1(j) = a1(j) + 72.0*cos(thta(j))
      a2(j) = a2(j) + 72.0*sin(thta(j))
    enddo
    M = M + 72.0
  enddo
  do j=1,3
    a1(j) = a1(j)/M
    a2(j) = a2(j)/M
    thta_v(j) = atan2(-a2(j),-a1(j)) + pi
    tmpcom(j) = L(j)*thta_v(j)/(2*pi)
    if (tmpcom(j) > L(j)/2) then
      tmpcom(j) = tmpcom(j) - L(j)
    endif
  enddo

End Subroutine

Subroutine Read_Input(atoms_per_species)
  ! Subroutine to read an input file
  use parameters
  use constants
  implicit none
  integer :: i
  integer, dimension(10) :: atoms_per_species
  open(10,file="atoms.in", status='old')
  read(10,*)
  read(10,*) nspecies
  read(10,*)
  read(10,*) nframes
  read(10,*)
  natoms = 0
  do i=1,nspecies
    read(10,*) atoms_per_species(i)
    natoms = natoms + atoms_per_species(i)
  enddo
End Subroutine

Subroutine Read_Frame(r)
  use constants
  use parameters
  implicit none
  integer :: i,j,k
  real, dimension(natoms,3) :: r
  read(ufile,*)
  read(ufile,*)
  read(ufile,*)
  read(ufile,*)
  read(ufile,*)
  read(ufile,*)
  read(ufile,*) xlo, xhi
  read(ufile,*) ylo, yhi
  read(ufile,*) zlo, zhi
  read(ufile,*)
  L(1) = xhi-xlo
  L(2) = yhi-ylo
  L(3) = zhi-zlo
  do i=1,natoms
    read(ufile,*) ctmp, ctmp, r(i,1), r(i,2), r(i,3)
  enddo 
End Subroutine
  
