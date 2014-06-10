!====================================================================
!   initial conditions
!   for spherical blast wave
!====================================================================
subroutine initflow(itprint)
  use parameters
  use globals
  use user_mod


  implicit none
  integer , intent(inout) :: itprint
  integer ::  i,j,k,unitin,err
  character (len=128) file1
  !   for the random number generator
  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  character (len=10) :: system_time
  real :: rtime,ran
  real :: x, y, z
  real :: rads, dens, velx, vely, velz, cpi
  !
  !
  call random_seed(size=rand_size)
  allocate(rand_seed(1:rand_size))
  call date_and_time(time=system_time)
  read(system_time,*) rtime
  rand_seed=int(rtime*1000.)
#ifdef MPIP
  rand_seed=rand_seed*rank
#endif
  call random_seed(put=rand_seed)    
  !
  !--------------------------------------------------------------------
  !
  if (iwarm.eq.0) then
    !
    call initial_conditions(u,0.)
    ! 
  else
     !
     !   read from previous (.bin) output
#ifdef MPIP
     write(file1,'(a,i3.3,a,i3.3,a)')  &
          trim(outputpath)//'BIN/points',rank,'.',itprint,'.bin'
     unitin=rank+10
#else
     write(file1,'(a,i3.3,a)')         &
          trim(outputpath)//'BIN/points',itprint,'.bin'
     unitin=10
#endif
     open(unit=unitin,file=file1,status='unknown',form='unformatted', &
          convert='LITTLE_ENDIAN')
     !
     read(unitin) u(:,:,:,:)
     close(unitin)
     !
     print'(i3,a,a)',rank,'read',trim(file1)
     itprint=itprint+1
     !
     !
#ifdef MPIP
     call mpi_barrier(mpi_comm_world,err)
#endif
     !
  end if
  !---------------------------------------------------------------------
  return
end subroutine initflow
!=======================================================================
