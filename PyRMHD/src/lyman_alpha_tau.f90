!=======================================================================
!   Computes the opacity of the Lyman Alpha line in nvmap velocity
!   chanels from vmin to vmax
!=======================================================================
program lyman_alpha_tau
  use parameters
  use globals, only : u, dx, dy, dz, coords, rank, left, right   &
                     , top, bottom, out, in, rank, comm2d
  implicit none
  character (len=128) :: filepath
  real, allocatable :: ph(:,:,:)
  integer :: err
  integer :: itprint
  ! 
  real, parameter :: theta_x = 3.33 *pi/180.
  real, parameter :: theta_y = 0.00 *pi/180.
  real, parameter :: theta_z = 0.00 *pi/180.
  !   map and its dimensions
  integer, parameter :: nxmap=300, nymap=300,nvmap=250
  real :: dxT, dyT, vmin,vmax
  real :: map(nxmap, nymap,nvmap), map1(nxmap, nymap,nvmap)
  !--------------------------------------------------------------------
  ! initializes program
  call init_LA()
  allocate( ph(nx,ny,nz) )
  !  minumim and maximum of the velocity channels
  vmin=-320.e5
  vmax= 200.e5
  !  data pixes size, relative to the simulation
  dxT=dx
  dyT=dy
  !
  ! chose output (fix later to input form screen)
  filepath='./0.3G-G1.01-RD-300/'

  loop_over_outputs : do itprint=0,38
    
    !  read ph and u from file
    call readdata(u,ph,itprint,filepath)
    !
    !  resets map
    map(:,:,:)=0.
    map1(:,:,:)=0.
    !
    if (rank == master) then
       print'(a)', 'Calculating projection with angles of rotaton'
       print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, ' &
                                    ,theta_y*180./pi,'° around Y, '&
                                    ,theta_z*180./pi,'° around Z, '
    end if
    !
    !  add info to the map
    call fillmap(nxmap,nymap,nvmap,vmin,vmax,u,map,dxT,dyT, theta_x, theta_y, theta_z)
    !  sum all the partial sums
    call mpi_reduce(map,map1,nxmap*nymap*nvmap, mpi_real_kind, mpi_sum, master, comm2d, err)
    
    !  write result
    if (rank == master) then 
      call write_LA(itprint,filepath,nxmap,nymap,nvmap,map1)
    end if
    
  end do loop_over_outputs

  if (rank == master) print*, 'my work here is done, have a  nice day'
#ifdef MPIP
  call mpi_finalize(err)
#endif
  !
  stop
end program lyman_alpha_tau
!=======================================================================

!=======================================================================
!  Initializes MPI, data arrays, etc
!=======================================================================
subroutine init_LA()
use parameters
use globals, only : u, dx, dy, dz, coords, rank, left, right   &
                     , top, bottom, out, in, rank, comm2d
implicit none
  integer :: nps, err
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period

  !initializes MPI
#ifdef MPIP
#ifdef PERIODX
  logical, parameter :: perx=.true.
#else
  logical, parameter :: perx=.false.
#endif
#ifdef PERIODY
  logical, parameter :: pery=.true.
#else
  logical, parameter :: pery=.false.
#endif
#ifdef PERIODZ
  logical, parameter :: perz=.true.
#else
  logical, parameter :: perz=.false.
#endif
  period(0)=perx
  period(1)=pery
  period(2)=perz
  dims(0)  =mpicol
  dims(1)  =mpirow
  dims(2)  =mpirowz
  !
  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)
  if (nps.ne.np) then
     print*, 'processor number (',nps,') is not equal to pre-defined number (',np,')'
     call mpi_finalize(err) 
     stop
  endif
#else
  rank=0
  coords(:)=0
#endif
  if(rank.eq.master) then
     print '(a)' ,"*******************************************"
     print '(a)' ," _                      _                 *"
     print '(a)' ,"| |__  _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ,"| '_ \| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| | | | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ,"|_| |_|\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ,"                                          *"
  endif
#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', np , ' processors    *'
     print '(a)' ,'*******************************************'
     print '(a)', 'Calculating Lyman Alpha Tau'
  end if
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, 1            &
       , comm2d, err)
  call mpi_comm_rank(comm2d, rank, err)
  call mpi_cart_coords(comm2d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                              &
       ,' ready w/coords',coords(0),coords(1),coords(2)   
  call mpi_cart_shift(comm2d, 0, 1, left  , right, err)
  call mpi_cart_shift(comm2d, 1, 1, bottom, top  , err)
  call mpi_cart_shift(comm2d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)   
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
  print '(a)', 'Calculating Lyman Alpha Tau'
#endif

!   grid spacing
  dx=xmax/nxtot
  dy=ymax/nytot
  dz=zmax/nztot

  !   allocate big arrays in memory
  allocate( u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

end subroutine init_LA
!=======================================================================

!=======================================================================
subroutine readdata(u,ph,itprint,filepath)
  use parameters, only : nx, ny, nz, np, &
                  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals, only : rank, comm2d
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax), &
                       ph(nx,ny,nz)
  integer, intent(in) :: itprint
  character (len=128), intent(in) :: filepath
  integer :: unitin, ip, err
  character (len=128) file1
  
  take_turns : do ip=0,np-1
    if (rank == ip) then
  !  Read u
#ifdef MPIP
        write(file1,'(a,i3.3,a,i3.3,a)')  &
             trim(filepath)//'BIN/points',rank,'.',itprint,'.bin'
        unitin=rank+10
#else
         write(file1,'(a,i3.3,a)')         &
              trim(filepath)//'BIN/points',itprint,'.bin'
         unitin=10
#endif
         open(unit=unitin,file=file1,status='unknown',form='unformatted', &
              convert='LITTLE_ENDIAN')
         !
         read(unitin) u(:,:,:,:)
         close(unitin)
         !
         print'(i3,a,a)',rank,' read file:',trim(file1)

      ! read ph (this assumes that MPI is enabled)
      write(file1,'(a,i3.3,a,i3.3,a)')  &
          trim(filepath)//'BIN/ph-',rank,'.',itprint,'.bin'
      unitin=10+rank

      open(unit=unitin,file=file1,status='unknown',form='unformatted', &
           convert='LITTLE_ENDIAN')
      read (unitin) ph(:,:,:)
      close(unitin)
      print'(i3,a,a)',rank," read file:",trim(file1)

    end if
    call mpi_barrier(comm2d,err)
  end do take_turns

end subroutine readdata
!=======================================================================

!=======================================================================
  !   gets the position and spherical radius calculated with respect to
  !   the center of the grid
  subroutine getXYZ(i,j,k,x,y,z)
    use globals,    only : dx, dy, dz, coords
    use parameters, only : nx, ny, nz, nxtot, nytot, nztot
    implicit none
    integer, intent(in)  :: i, j, k
    real,    intent(out) :: x, y, z
 
    x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
    y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
    z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
        
  end subroutine getXYZ
!=======================================================================

!=======================================================================
subroutine rotation_x(theta,x,y,z,xn,yn,zn)
   ! rotation around the x axis by an angle theta
   implicit none
   real, intent(in ) :: theta, x, y, z
   reaL, intent(out) :: xn, yn, zn
   xn =   x
   yn =   y*cos(theta) - z*sin(theta)
   zn =   y*sin(theta) + z*cos(theta)
 end subroutine rotation_x
!=======================================================================

!=======================================================================
 subroutine rotation_y(theta,x,y,z,xn,yn,zn)
   ! rotation around the y axis by an angle theta
   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) + z*sin(theta)
   yn =   y
   zn = - x*sin(theta) + z*cos(theta)
 end subroutine rotation_y
!=======================================================================

!=======================================================================
 subroutine rotation_z(theta,x,y,z,xn,yn,zn)
   ! rotation around the z axis by an angle theta
   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) - y*sin(theta)
   yn =   x*sin(theta) + y*cos(theta)
   zn =   z
 end subroutine rotation_z
!=======================================================================

!=======================================================================
subroutine fillmap(nxmap,nymap,nvmap,vmin,vmax,u,map,dxT,dyT,theta_x,theta_y,theta_z)
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc2, rsc, clight,nztot
  use globals, only : rank, dz
  implicit none
  integer, intent(in) :: nxmap,nymap,nvmap
  real, intent(in) :: vmin, vmax 
  real, intent(in) :: u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real , intent(in) :: dxT, dyT, theta_x, theta_y, theta_z
  real, intent(out) :: map(nxmap,nymap,nvmap)
  integer :: i,j,k, iobs, jobs
  real :: x,y,z,xn,yn,zn, xobs, yobs, vx, vy, vz,vxn, vyn, vzn, velsc, z0
  real :: T, prim(neq), profile(nvmap)
  real, parameter :: sigmaLA = 0.01105, lambdaLA=1.215668e-5 !(c/nu0=lambda)
  velsc=sqrt(vsc2)

  do i=1,nx
     do j=1,ny
        do k=1, nz
          !  obtain original position
          call getXYZ(i,j,k, x,y,z)
          !only do to the side facing the observer
          if(z < 0) then
            !  do the rotation of the coordinates
           call rotation_x(theta_x,x,y,z,xn,yn,zn)
            call rotation_y(theta_y,xn,yn,zn,x,y,z)
            call rotation_z(theta_z,x,y,z,xn,yn,zn)
            ! This is the position on the target (centered)
            ! Integration is along Z
            iobs=xn/dxT + nxmap/2
            jobs=yn/dyT + nymap/2
            !----------------------------------------
            !  get the velocity in cm/s
            call uprim(prim, u(:,i,j,k),T)
            vx=prim(2)*velsc
            vy=prim(3)*velsc
            vz=prim(4)*velsc
            !  obtain the LOS velocity
            call rotation_x(theta_x,vx,vy,vz,vxn,vyn,vzn)
            call rotation_y(theta_y,vxn,vyn,vzn,vx,vy,vz)
            call rotation_z(theta_z,vx,vy,vz,vxn,vyn,vzn)
            !----------------------------------------
            !  calculate the line profile function
            call phigauss(T,vzn,vmin,vmax,nvmap,profile) 
            !  make sure the result lies in the map bounds
            if( (iobs >=1    ).and.(jobs >=1    ).and. &
                (iobs <=nxmap).and.(jobs <=nymap) ) then
              !if ((T < 1e5).and.(prim(7)<0)) then
                map(iobs,jobs,:)= map(iobs,jobs,:) + &
                                  dz*rsc*prim(neqdin+1)*sigmaLA*lambdaLA*profile(:)
              !end if
            end if
          end if
      end do
    end do
  end do

end subroutine fillmap
!=======================================================================

!=======================================================================
subroutine  write_LA(itprint,filepath,nxmap,nymap,nvmap,map)
  implicit none
  integer, intent(in) :: nxmap, nymap,nvmap,itprint
  character (len=128), intent(in) :: filepath
  real, intent(in) :: map(nxmap,nymap,nvmap)
  character (len=128) file1
  integer ::  i, j, k, unitout, ip

  write(file1,'(a,i3.3,a)')  trim(filepath)//'BIN/LA_tau-',itprint,'.bin'
  unitout=11
  open(unit=unitout,file=file1,status='unknown',form='unformatted', &
       convert='LITTLE_ENDIAN')
  
  write (unitout) map(:,:,:)
  close(unitout)
  print'(a,a)'," wrote file:",trim(file1)
  
end subroutine write_LA
!=======================================================================

!=======================================================================
!  This routine computes a gaussian line profile
!=======================================================================
subroutine phigauss(T,vzn,vmin,vmax,nvmap,profile) 
  use parameters, only: amh, pi, kB, clight
  implicit none
  real, intent(in) :: T, vzn, vmin, vmax
  integer, intent(in) :: nvmap
  real, intent(out) :: profile(nvmap)
  integer :: i
  real :: coef, dv, vr
  
  profile(:)=0.
  dv=(vmax-vmin)/float(nvmap)
  
  coef=amh/(2.*kB*T)
  
  do i=1,nvmap
     vr=(float(i)-0.5)*dv+vmin
     profile(i)=sqrt(coef/pi)*exp(-coef*((vr-vzn)**2) )
  end do
  
end subroutine phigauss
!=======================================================================

