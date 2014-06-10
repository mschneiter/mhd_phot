!======================================================================
!   boundary conditions on up (two ghost cells)
!======================================================================
subroutine boundaryII(time,dt)
  use parameters
  use globals
#ifdef OTHERB
!  use jet
  use star
#endif
  implicit none
  real, intent(in) :: time,dt
  integer :: i, j, k, l,ns
  integer, parameter :: nxmg=nx-nghost+1 ,nxp=nx+1
  integer, parameter :: nymg=ny-nghost+1, nyp=ny+1
  integer, parameter :: nzmg=nz-nghost+1, nzp=nz+1
  real :: ekin,x,y,z,r,rho!r1, r2, rx, ry, rz, vjet, xp, yp, zp
  real, allocatable::xc(:),yc(:),zc(:) !star positions
  !real :: omegat1, omegat2
  !real :: coso, coso2, sino, sino2, cosasino,cosacoso,sinasino,sinacoso
#ifdef MPIP
  integer:: status(MPI_STATUS_SIZE), err
  real, dimension(neq,nghost,nymin:nymax,nzmin:nzmax)::sendr,recvr,sendl,recvl
  real, dimension(neq,nxmin:nxmax,nghost,nzmin:nzmax)::sendt,recvt,sendb,recvb
  real, dimension(neq,nxmin:nxmax,nymin:nymax,nghost)::sendi,recvi,sendo,recvo
  integer, parameter :: bxsize=neq*nghost*(nymax-nymin+1)*(nzmax-nzmin+1)
  integer, parameter :: bysize=neq*(nxmax-nxmin+1)*nghost*(nzmax-nzmin+1)
  integer, parameter :: bzsize=neq*(nxmax-nxmin+1)*(nymax-nymin+1)*nghost
  !
#endif
  !
#ifdef MPIP
     !
     !   Exchange boundaries between processors
     !   -------------------------------------------------------------
     !
     !   boundaries to processors to the right, left, top, and bottom
     sendr(:,1:nghost,:,:)=up(:,nxmg:nx ,:,:)
     sendl(:,1:nghost,:,:)=up(:,1:nghost,:,:)
     sendt(:,:,1:nghost,:)=up(:,:,nymg:ny ,:)
     sendb(:,:,1:nghost,:)=up(:,:,1:nghost,:)
     sendi(:,:,:,1:nghost)=up(:,:,:,nzmg:nz )
     sendo(:,:,:,1:nghost)=up(:,:,:,1:nghost)
     !
     call mpi_sendrecv(sendr, bxsize, mpi_real_kind, right  ,0,          &
                       recvl, bxsize, mpi_real_kind, left   ,0,          &
                       comm2d, status , err)
     !
     call mpi_sendrecv(sendt, bysize, mpi_real_kind, top    ,0,          &
                       recvb, bysize, mpi_real_kind, bottom ,0,          &
                       comm2d, status , err)
     !
     call mpi_sendrecv(sendi, bzsize, mpi_real_kind, in     ,0,          &
                       recvo, bzsize, mpi_real_kind, out    ,0,          &
                       comm2d, status , err)
     !---
     call mpi_sendrecv(sendl, bxsize, mpi_real_kind, left  , 0,          &
                       recvr, bxsize, mpi_real_kind, right , 0,          &
                       comm2d, status , err)
     !
     call mpi_sendrecv(sendb, bysize, mpi_real_kind, bottom, 0,          &
                       recvt, bysize, mpi_real_kind, top   , 0,          &
                       comm2d, status , err)
     !
     call mpi_sendrecv(sendo, bzsize, mpi_real_kind, out   , 0,          &
                       recvi, bzsize, mpi_real_kind, in    , 0,          &
                       comm2d, status , err)
     !
     if (left  .ne. -1) up(:,nxmin:0  ,:,:)=recvl(:,1:nghost,:,:)
     if (right .ne. -1) up(:,nxp:nxmax,:,:)=recvr(:,1:nghost,:,:)
     if (bottom.ne. -1) up(:,:,nymin:0  ,:)=recvb(:,:,1:nghost,:)
     if (top   .ne. -1) up(:,:,nyp:nymax,:)=recvt(:,:,1:nghost,:)
     if (out   .ne. -1) up(:,:,:,nzmin:0  )=recvo(:,:,:,1:nghost)
     if (in    .ne. -1) up(:,:,:,nzp:nzmax)=recvi(:,:,:,1:nghost)
     !
#else
     !
     !******************************************************************
     !   periodic BCs
#ifdef PERIODX
     !-----------------------------------
     !   Left BC
     if (coords(0).eq.0) then
        up(:,nxmin:0,:,:)=up(:,nxmg:nx,:,:)
     endif
     !   Right BC
     if (coords(0).eq.mpicol-1) then
        up(:,nxp:nxmax,:,:)=up(:,1:nghost,:,:)
     endif
     !-----------------------------------
#endif
#ifdef PERIODY
     !-----------------------------------
     !   bottom BC
     if (coords(1).eq.0) then
        up(:,:,nymin:0,:)= up(:,:,nymg:ny,:)
     endif
     !   top BC
     if (coords(1).eq.mpirow-1) then
        up(:,:,nyp:nymax,:)= up(:,:,1:nghost,:)
     endif
     !-----------------------------------
#endif
#ifdef PERIODZ
     !-----------------------------------
     !   out BC
     if (coords(2).eq.0) then
        up(:,:,:,nzmin:0)= up(:,:,:,nzmg:nz)
     endif
     !   in BC
     if (coords(2).eq.mpirowz-1) then
        up(:,:,:,nzp:nzmax)= up(:,:,:,1:nghost)
     endif
     !-----------------------------------
#endif
     !******************************************************************
     !
#endif   !MPIP
     !
     !******************************************************************
     !   Reflecting BCs
#ifdef REFXL
     if (coords(0).eq.0) then
        j=nghost
        do i=nxmin,0
           up(1  ,i,:,:) = up(1  ,j,:,:)
           up(2  ,i,:,:) =-up(2  ,j,:,:)
           up(3:neq,i,:,:) = up(3:neq,j,:,:)
           j=j-1
        enddo
     endif
#endif
     !-----------------------------------
#ifdef REFXR
     if (coords(0).eq.mpicol-1) then
        j=nx
        do i=nxp,nxmax
           up(1  ,i,:,:) = up(1  ,j,:,:)
           up(2  ,i,:,:) =-up(2  ,j,:,:)
           up(3:neq,i,:,:) = up(3:neq,j,:,:)
           j=j-1
        enddo
     endif
#endif
     !-----------------------------------
#ifdef REFYB
     if (coords(1).eq.0) then
        j=nghost
        do i=nymin,0
           up(1:2,:,i,:) = up(1:2,:,j,:)
           up(3  ,:,i,:) =-up(3  ,:,j,:)
           up(4:neq,:,i,:) = up(4:neq,:,j,:)
           j=j-1
        enddo
     endif
#endif
!-----------------------------------
#ifdef REFYT
     if (coords(1).eq.(mpirow-1)) then
        j=ny
        do i=nyp,nymax
           up(1:2,:,i,:) = up(1:2,:,j,:)
           up(3  ,:,i,:) =-up(3  ,:,j,:)
           up(4:neq,:,i,:) = up(4:neq,:,j,:)
           j=j-1
        enddo
     endif
#endif
!-----------------------------------
#ifdef REFZO
     if (coords(2).eq.0) then
        j=nghost
        do i=nzmin,0
           up(1:3,:,:,i) = up(1:3,:,:,j)
           up(4  ,:,:,i) =-up(4  ,:,:,j)
           up(5:neq  ,:,:,i) = up(5:neq  ,:,:,j)
           j=j-1
        enddo
     endif
#endif
!-----------------------------------
#ifdef REFZI
     if (coords(2).eq.mpirowz-1) then
        j=nz
        do i=nzp,nzmax
           up(1:3,:,:,i) = up(1:3,:,:,j)
           up(4  ,:,:,i) =-up(4  ,:,:,j)
           up(5:neq  ,:,:,i) = up(5:neq  ,:,:,j)
           j=j-1
        enddo
     endif
#endif
     !******************************************************************
     !
     !******************************************************************
     !   outflow BCs
     !-----------------------------------
     !   left
#ifdef OUTFXL
     if (coords(0).eq.0) then
        do i=nxmin,0
           up(:,i,:,:)=up(:,1,:,:)
        enddo
     endif
#endif
     !-----------------------------------
     !   right
#ifdef OUTFXR
     if (coords(0).eq.mpicol-1) then
        do i=nxp,nxmax
           up(:,i,:,:)=up(:,nx,:,:)
        enddo
     endif
#endif
     !-----------------------------------
     !   bottom
#ifdef OUTFYB
     if (coords(1).eq.0) then
        do i=nymin,0
           up(:,:,i,:)=up(:,:,1,:)
        enddo
     endif
#endif
     !-----------------------------------
     !   top
#ifdef OUTFYT
     if (coords(1).eq.mpirow-1) then
        do i=nyp,nymax
           up(:,:,i,:)=up(:,:,ny,:)
        enddo
     endif
#endif
     !-----------------------------------
     !   out
#ifdef OUTFZO
     if (coords(2).eq.0) then
        do i=nzmin,0
           up(:,:,:,i)=up(:,:,:,1)
        enddo
     endif
#endif
     !-----------------------------------
     !   in
#ifdef OUTFZI
     if (coords(2).eq.mpirowz-1) then
        do i=nzp,nzmax
           up(:,:,:,i)=up(:,:,:,nz)
        enddo
     endif
#endif
     !
     !******************************************************************
     !
     !******************************************************************
     !   other type of bounadries  <e.g. winds jets outflows>
#ifdef OTHERB
     !
     !   impose star/planet system (parameteres in user_mod.f90)
     call impose_user_bc(up,time)
     !
     !--------------------------------------------------------------------------------
#endif
  !
  return
end subroutine boundaryII
!//////////////////////////////////////////////////////////////////////
