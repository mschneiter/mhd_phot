!======================================================================
!   boundary conditions on u with only one ghost cell
!======================================================================
subroutine boundaryI(time,dt)
  use parameters
  use globals
#ifdef OTHERB
!  use jet
  use star
#endif
  implicit none
  real, intent(in) :: time,dt
  integer :: i,j,k,l
  integer, parameter :: nxm1=nx-1 ,nxp1=nx+1
  integer, parameter :: nym1=ny-1, nyp1=ny+1
  integer, parameter :: nzm1=nz-1, nzp1=nz+1
  real :: ekin,x,y,z,r,rho
  !real :: omegat1, omegat2
  !real :: coso, coso2, sino, sino2, cosasino,cosacoso,sinasino,sinacoso
  real, allocatable :: xc(:), yc(:), zc(:) !star positions
#ifdef MPIP
  integer:: status(MPI_STATUS_SIZE), err
  real, dimension(neq,1,0:nyp1,0:nzp1)::sendr,recvr,sendl,recvl
  real, dimension(neq,0:nxp1,1,0:nzp1)::sendt,recvt,sendb,recvb
  real, dimension(neq,0:nxp1,0:nyp1,1)::sendi,recvi,sendo,recvo
  integer, parameter :: bxsize=neq*(ny+2)*(nz+2)
  integer, parameter :: bysize=neq*(nx+2)*(nz+2)
  integer, parameter :: bzsize=neq*(nx+2)*(ny+2)
  !
#endif
  !
#ifdef MPIP
     !
     !   Exchange boundaries between processors
     !   -------------------------------------------------------------
     !
     !   boundaries to procs: right, left, top, bottom, in and out
     sendr(:,1,:,:)=u(:,nx    ,0:nyp1,0:nzp1)
     sendl(:,1,:,:)=u(:,1     ,0:nyp1,0:nzp1)
     sendt(:,:,1,:)=u(:,0:nxp1,ny    ,0:nzp1)
     sendb(:,:,1,:)=u(:,0:nxp1,1     ,0:nzp1)
     sendi(:,:,:,1)=u(:,0:nxp1,0:nyp1,nz    )
     sendo(:,:,:,1)=u(:,0:nxp1,0:nyp1,1     )
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
     if (left  .ne. -1) u(:,0     ,0:nyp1,0:nzp1)=recvl(:,1,:,:)
     if (right .ne. -1) u(:,nxp1  ,0:nyp1,0:nzp1)=recvr(:,1,:,:)
     if (bottom.ne. -1) u(:,0:nxp1,0     ,0:nzp1)=recvb(:,:,1,:)
     if (top   .ne. -1) u(:,0:nxp1,nyp1  ,0:nzp1)=recvt(:,:,1,:)
     if (out   .ne. -1) u(:,0:nxp1,0:nyp1,0     )=recvo(:,:,:,1)
     if (in    .ne. -1) u(:,0:nxp1,0:nyp1,nzp1  )=recvi(:,:,:,1)
     !
#else
     !
     !******************************************************************
     !   periodic BCs
#ifdef PERIODX
     !-----------------------------------
     !   Left BC
     if (coords(0).eq.0) then
        u(:,0,:,:)=u(:,nx,:,:)
     endif
     !   Right BC
     if (coords(0).eq.mpicol-1) then
        u(:,nxp1,:,:)=u(:,1,:,:)
     endif
     !-----------------------------------
#endif
#ifdef PERIODY
     !-----------------------------------
     !   bottom BC
     if (coords(1).eq.0) then
        u(:,:,0,:)= u(:,:,ny,:)
     endif
     !   top BC
     if (coords(1).eq.mpirow-1) then
        u(:,:,nyp1,:)= u(:,:,1,:)
     endif
     !-----------------------------------
#endif
#ifdef PERIODZ
     !-----------------------------------
     !   out BC
     if (coords(2).eq.0) then
        u(:,:,:,0)= u(:,:,:,nz)
     endif
     !   in BC
     if (coords(2).eq.mpirowz-1) then
        u(:,:,:,nzp1)= u(:,:,:,1)
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
        u(1       ,0,0:nyp1,0:nzp1) = u(1       ,1,0:nyp1,0:nzp1)
        u(2       ,0,0:nyp1,0:nzp1) =-u(2       ,1,0:nyp1,0:nzp1)
        u(3:neq,0,0:nyp1,0:nzp1) = u(3:neq,1,0:nyp1,0:nzp1)
     endif
#endif
     !-----------------------------------
#ifdef REFXR
     if (coords(0).eq.(mpicol-1)) then
        u(1       ,nxp1,0:nyp1,0:nzp1) = u(1       ,nx,0:nyp1,0:nzp1)
        u(2       ,nxp1,0:nyp1,0:nzp1) =-u(2       ,nx,0:nyp1,0:nzp1)
        u(3:neq,nxp1,0:nyp1,0:nzp1) = u(3:neq,nx,0:nyp1,0:nzp1)
     endif
#endif
     !-----------------------------------
#ifdef REFYB
     if (coords(1).eq.0) then
        u(1:2     ,0:nxp1,0,0:nzp1) = u(1:2     ,0:nxp1,1,0:nzp1)
        u(3       ,0:nxp1,0,0:nzp1) =-u(3       ,0:nxp1,1,0:nzp1)
        u(4:neq,0:nxp1,0,0:nzp1) = u(4:neq,0:nxp1,1,0:nzp1)
     endif
#endif
!-----------------------------------
#ifdef REFYT
     if (coords(1).eq.(mpirow-1)) then
        u(1:2     ,0:nxp1,nyp1,0:nzp1) = u(1:2     ,0:nxp1,ny,0:nzp1)
        u(3       ,0:nxp1,nyp1,0:nzp1) =-u(3       ,0:nxp1,ny,0:nzp1)
        u(4:neq,0:nxp1,nyp1,0:nzp1) = u(4:neq,0:nxp1,ny,0:nzp1)
     endif
#endif
!-----------------------------------
#ifdef REFZO
     if (coords(2).eq.0) then
           u(1:3     ,0:nxp1,0:nyp1,0) = u(1:3     ,0:nxp1,0:nyp1,1)
           u(4       ,0:nxp1,0:nyp1,0) =-u(4       ,0:nxp1,0:nyp1,1)
           u(5:neq,0:nxp1,0:nyp1,0) = u(5:neq,0:nxp1,0:nyp1,1)
     endif
#endif
!-----------------------------------
#ifdef REFZI
     if (coords(2).eq.mpirowz-1) then
           u(1:3     ,0:nxp1,0:nyp1,nzp1) = u(1:3     ,0:nxp1,0:nyp1,nz)
           u(4       ,0:nxp1,0:nyp1,nzp1) =-u(4       ,0:nxp1,0:nyp1,nz)
           u(5:neq,0:nxp1,0:nyp1,nzp1) = u(5:neq,0:nxp1,0:nyp1,nz)
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
        u(:,0,   0:nyp1,0:nzp1)=u(:,1 ,0:nyp1,0:nzp1)
     endif
#endif
     !-----------------------------------
     !   right
#ifdef OUTFXR
     if (coords(0).eq.mpicol-1) then
        u(:,nxp1,0:nyp1,0:nzp1)=u(:,nx,0:nyp1,0:nzp1)
     endif
#endif
     !-----------------------------------
     !   bottom
#ifdef OUTFYB
     if (coords(1).eq.0) then
        u(:,0:nxp1,0   ,0:nzp1)=u(:,0:nxp1,1 ,0:nzp1)
     endif
#endif
     !-----------------------------------
     !   top
#ifdef OUTFYT
     if (coords(1).eq.mpirow-1) then
        u(:,0:nxp1,nyp1,0:nzp1)=u(:,0:nxp1,ny,0:nzp1)
     endif
#endif
     !-----------------------------------
     !   out
#ifdef OUTFZO
     if (coords(2).eq.0) then
           u(:,0:nxp1,0:nyp1,0   )=u(:,0:nxp1,0:nyp1,1 )
     endif
#endif
     !-----------------------------------
     !   in
#ifdef OUTFZI
     if (coords(2).eq.mpirowz-1) then
           u(:,0:nxp1,0:nyp1,nzp1)=u(:,0:nxp1,0:nyp1,nz)
     endif
#endif
     !
     !******************************************************************
     !
     !******************************************************************
     !   other type of boundaries
#ifdef OTHERB
     !
     ! imposte star/planet system (parameters in user_mod.f90)
    call impose_user_bc(u,time)
    !
    !--------------------------------------------------------------------------------
    !
#endif
  !
  return
end subroutine boundaryI
!//////////////////////////////////////////////////////////////////////
