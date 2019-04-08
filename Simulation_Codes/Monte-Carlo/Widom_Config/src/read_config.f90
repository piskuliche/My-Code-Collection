! Subroutine for reading data

Subroutine read_data(config,skip_config,noconfigs,nmol)

      implicit none
      integer :: i, j, k
      integer :: skip_config
      integer :: nmol
      character(len=10) config
      real, dimension(:,:,:) :: ALLOCATABLE :: r
        
      ALLOCATE(r(noconfigs,nmol,3))
      
      open(10,file=trim(config),status='old')

      ! Read in configuration
      read(10,*) natms
      read(10,*)
      do i = 1, nmol
        do k = 1, atms_per_mol
            read(10,*) 
        enddo
      enddo

