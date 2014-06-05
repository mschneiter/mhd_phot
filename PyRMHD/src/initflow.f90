!====================================================================
!   initial conditions
!   for spherical blast wave
!====================================================================
subroutine initflow(itprint)
  use parameters
  use globals
#ifdef OTHERB
   use star
  ! use jet
#endif
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
     do i=nxmin,nxmax
        do j=nymin,nymax
           do k=nzmin,nzmax
             ! Position measured from the centre of the grid (star)              
             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
             z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

             ! Distance from the center of the box
             rads=sqrt(x**2+y**2+z**2)
             if(rads == 0.) rads=dx*0.01

             !  This is the solution of the wind of the star
             !  filling the whole domain
             cpi=b0*(rsw/(rads+1.e-30))**3/(2.*(rads+1.e-30)**2)
             
             VelX=VSW*X/RADS
             VelY=VSW*Y/RADS
             VelZ=VSW*Z/RADS
             DENS=DSW*RSW**2/RADS**2
             !   total density and momenta
             u(1,i,j,k) = dens
             u(2,i,j,k) = dens*velx
             u(3,i,j,k) = dens*vely
             u(4,i,j,k) = dens*velz
             ! B FIELD
#ifdef PMHD
             u(6,i,j,k) =  3.*y*x*cpi
             u(7,i,j,k) = (3.*y**2-rads**2)*cpi
             u(8,i,j,k) =  3.*y*z*cpi
#elif MHD
             u(6,i,j,k) =  3.*y*x*cpi
             u(7,i,j,k) = (3.*y**2-rads**2)*cpi
             u(8,i,j,k) =  3.*y*z*cpi
#endif 
             !   energy (5)
             u(5,i,j,k)=0.5*dens*vsw**2 &
                  + cv*dens*2.*Tsw

             !  density of neutrals (6) ahora (9)
             u(neqdin+1,i,j,k)= 0.*dens
             !   passive scalar (7) ahora (10)
             u(neqdin+2,i,j,k)= dens
             !   
           end do
        end do
     end do
     !
     !---------------------------------------------------------------------

     !
     !
     !*****************************************************************
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
