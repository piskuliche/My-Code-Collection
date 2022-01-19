! Copyright 2022, Zeke A. Piskulich, Boston University
! This code takes gREM trajectories and sorts them into walker files (walker-XXX.lammpsdump)
!
! There is one required input file: grem.input found in the current directory.
! This code outputs the walkers to a new directory.



module constants
    implicit none
    real, parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: numfile
  integer :: fstart, fstop, nwindows,natoms, nlipids
  integer :: nlog, dumpfreq
  real, dimension(3) :: L
  common fstart, fstop, nwindows,natoms,nlog,dumpfreq,nlipids
  common L
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
  integer :: outvalue
  integer, dimension(200) :: sortdata, infiles, outfiles

  character(len=20) :: str,ctmp

  infiles=0
  outfiles=0

  write(*,*) "Begin Program"
  call read_input()
  write(*,*) fstart,fstop,nwindows,nlog


  ! Open Output Files
  do k=1,nwindows
    outfiles(k) = 300 + k
    open(outfiles(k), file="all_dumps/walker-"//trim(str(k-1))//".lammpsdump")
  enddo !k
  

  do i=fstart, fstop

    ! Open log and read first 3 lines
    open(12,file="log/log.lammps-"//trim(str(i)), status='old')
    read(12,*)
    read(12,*)
    read(12,*)
    
    ! Open window files to read
    do k=1,nwindows
      infiles(k) = 100 + k
      open(infiles(k), file=trim(str(k-1))//"/dump-"//trim(str(i))//".dat",status='old')
    enddo

    ! loop over data in the log
    do j=1,nlog
      call readwalkers(sortdata)
      if (mod(j-1,dumpfreq) == 0) then
        do k=1,nwindows
          outvalue = outfiles(sortdata(k)+1)
          call readframe(infiles(k),outvalue)
        enddo !k
      endif
    enddo !j
    close(12)

    ! Close the Windows
    do k=1, nwindows
      close(infiles(k))
    enddo!k
  enddo !i

  ! Close Output files
  do k=1,nwindows
    close(outfiles(k))
  enddo
end program analyze_grem

subroutine readframe(infile, outfile)
  use constants
  use parameters
  implicit none
  integer :: i,j,k
  integer :: infile, outfile

  character(len=200) :: ctmp
  ! header
  do i=1,9
    read(infile,'(A)') ctmp
    write(outfile,*) trim(ctmp)
  enddo !i
  do i=1,nlipids
    do j=1,natoms
      read(infile,'(A)') ctmp
      write(outfile,*) trim(ctmp)
    enddo !j
  enddo !i
end subroutine readframe


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
  
  integer :: infile
  infile = 11
  open(infile, file="grem.input", status='old')
  read(infile,*)
  read(infile,*) nwindows
  read(infile,*)
  read(infile,*) fstart, fstop
  read(infile,*)
  read(infile,*) nlog, dumpfreq
  read(infile,*) 
  read(infile,*) natoms, nlipids
  close(infile)
end subroutine read_input
