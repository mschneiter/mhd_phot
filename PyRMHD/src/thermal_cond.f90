!=======================================================================
!  Module to add (isotropic) thermal conduction
!=======================================================================
#ifdef THERMAL_COND
module thermal_cond
  use globals
  use parameters
  implicit none
  real, parameter :: ph=0.4, nu=0.005, snu=sqrt(nu)
  real, allocatable :: Temp(:,:,:)
  real :: dtcond

contains
  !--------------------------------------------------------------------
  !  initializes array
  !--------------------------------------------------------------------
  subroutine init_thermal_cond
    implicit none

    allocate(Temp(nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

  end subroutine init_thermal_cond
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !   computes conduction timescale
  !--------------------------------------------------------------------
  subroutine dt_cond(dt)
    implicit none
    real, intent(out) :: dt
    real              :: dtp, cs, ddx
    integer :: i, j, k, err
    !
    dtp=1.e30
    ddx=min(dx ,dy)
    ddx=min(ddx,dz)
    !
    do i=1,nx
       do j=1,ny
          do k=1,nz
             !
             call csound(primit(5,i,j,k),primit(1,i,j,k),cs)
             !  spitzer timescale
             dtp=min( dtp, cv*primit(5,i,j,k)*Psc*(ddx*rsc)**2/(Ksp(Temp(i,j,k))*Temp(i,j,k) ) )
             !  saturated conduction timescale
             dtp=min( dtp, cv*ddx/(5*ph*cs)*tsc )
             !
          end do
       end do
    end do
    !

    dtp=dtp*0.25

#ifdef MPIP
    call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
    dt=dtp
#endif

  end subroutine dt_cond
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  progress bar (takes a number between 1 and tot)
  !--------------------------------------------------------------------
  subroutine progress(j,tot)
    implicit none
    integer(kind=4)::j,k
    integer(kind=4), intent(in) :: tot
    character(len=57)::bar="???% |                                                  |"
    open (unit=6, carriagecontrol='fortran')
    write(unit=bar(1:3),fmt="(i3)") 100*j/tot
    bar(7:56)="."
    do k=1, 50*j/tot
       bar(6+k:6+k)="="
    enddo
    ! print the progress bar.
    write(unit=6,fmt="(a1,a1,a57)") '+',char(13), bar
    return
  end subroutine progress
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !   Spitzer conductivity
  !--------------------------------------------------------------------
  real function KSp(T)
    implicit none
    real, intent(in) :: T
    real, parameter :: beta=6.e-7
    !
    Ksp= beta*T**(2.5)
    !
  end function KSp
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !   Heat flux, takes minimum of spitzer and saturated value
  !   The result is stored in the 5th component of global the
  !   F,G,H fluxes (in cgs) 
  !--------------------------------------------------------------------
  subroutine heatfluxes()
    use sound
    implicit none
    integer :: i, j, k
    real, parameter :: clight=3.E10, phi=0.3
    real:: cs, coef, dTx, dTy, dTz, meanT, meanP, meanDens, yhp
    !    
    do i=0,nx
       do j=0,ny
          do k=0,nz
             !
             yhp=1.-primit(neqdin+1,i,j,k)/primit(1,i,j,k)
             !
             !  get the flux in the X direction
             if (Temp(i,j,k) == Temp(i+1,j,k) ) then
                F(5,i,j,k)=0.
             else
                meanP   = 0.5*(primit(5,i,j,k)+primit(5,i+1,j,k))
                meanDens= 0.5*(primit(1,i,j,k)+primit(1,i+1,j,k))
                meanT   = 0.5*(    Temp(i,j,k)+    Temp(i+1,j,k))
                call csound(meanP,meanDens,cs)
                cs=min(cs*sqrt(vsc2),clight)

                dTx=(Temp(i+1,j,k)-Temp(i,j,k))/(dx*rsc)
                coef=min( Ksp(meanT) , 5.*ph*cs*meanP*Psc/abs(dTx) )

                F(5,i,j,k)=-coef*dTx*yhp

             end if
             !  get the flux in the Y direction
             if (Temp(i,j,k) == Temp(i,j+1,k) ) then
                G(5,i,j,k)=0.
             else
                meanP   = 0.5*(primit(5,i,j,k)+primit(5,i,j+1,k))
                meanDens= 0.5*(primit(1,i,j,k)+primit(1,i,j+1,k))
                meanT   = 0.5*(    Temp(i,j,k)+    Temp(i,j+1,k))
                call csound(meanP,meanDens,cs)
                cs=min(cs*sqrt(vsc2),clight)

                dTy=(Temp(i,j+1,k)-Temp(i,j,k))/(dy*rsc)   
                coef=min( Ksp(meanT) , 5.*ph*cs*meanP*Psc/abs(dTy) )
                G(5,i,j,k)=-coef*dTy*yhp

             end if
             !  get the flux in the Z direction
             if (Temp(i,j,k) == Temp(i,j,k+1) ) then
                H(5,i,j,k)=0.
             else
                meanP   = 0.5*(primit(5,i,j,k)+primit(5,i,j,k+1))
                meanDens= 0.5*(primit(1,i,j,k)+primit(1,i,j,k+1))
                meanT   = 0.5*(  Temp(i,j,k)  +    Temp(i,j,k+1))
                call csound(meanP,meanDens,cs)
                cs=min(cs*sqrt(vsc2),clight)

                dTz=(Temp(i,j,k+1)-Temp(i,j,k))/(dz*rsc)
                coef=min( Ksp(meanT) , 5.*ph*cs*meanP*Psc/abs(dTz) )

                H(5,i,j,k)=-coef*dTz*yhp

             end if

          end do
       end do
    end do
    !
  end subroutine heatfluxes
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  Exchanges one layer of boundaries, only the equation that
  !  corresponds to the energy
  !--------------------------------------------------------------------
  subroutine thermal_bounds()
    !use parameters
    !use globals
    implicit none
    integer, parameter :: nxp1=nx+1
    integer, parameter :: nyp1=ny+1
    integer, parameter :: nzp1=nz+1
#ifdef MPIP
    integer:: status(MPI_STATUS_SIZE), err
    real, dimension(1,1,0:nyp1,0:nzp1)::sendr,recvr,sendl,recvl
    real, dimension(1,0:nxp1,1,0:nzp1)::sendt,recvt,sendb,recvb
    real, dimension(1,0:nxp1,0:nyp1,1)::sendi,recvi,sendo,recvo
    integer, parameter :: bxsize=(ny+2)*(nz+2)
    integer, parameter :: bysize=(nx+2)*(nz+2)
    integer, parameter :: bzsize=(nx+2)*(ny+2)
    !
    !
    !   Exchange boundaries between processors
    !   -------------------------------------------------------------
    !
    !   boundaries to procs: right, left, top, bottom, in and out
    sendr(1,1,:,:)=u(5,nx    ,0:nyp1,0:nzp1)
    sendl(1,1,:,:)=u(5,1     ,0:nyp1,0:nzp1)
    sendt(1,:,1,:)=u(5,0:nxp1,ny    ,0:nzp1)
    sendb(1,:,1,:)=u(5,0:nxp1,1     ,0:nzp1)
    sendi(1,:,:,1)=u(5,0:nxp1,0:nyp1,nz    )
    sendo(1,:,:,1)=u(5,0:nxp1,0:nyp1,1     )
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
    if (left  .ne. -1) u(5,0     ,0:nyp1,0:nzp1)=recvl(1,1,:,:)
    if (right .ne. -1) u(5,nxp1  ,0:nyp1,0:nzp1)=recvr(1,1,:,:)
    if (bottom.ne. -1) u(5,0:nxp1,0     ,0:nzp1)=recvb(1,:,1,:)
    if (top   .ne. -1) u(5,0:nxp1,nyp1  ,0:nzp1)=recvt(1,:,1,:)
    if (out   .ne. -1) u(5,0:nxp1,0:nyp1,0     )=recvo(1,:,:,1)
    if (in    .ne. -1) u(5,0:nxp1,0:nyp1,nzp1  )=recvi(1,:,:,1)
    !
#else
    !
    !******************************************************************
    !   periodic BCs
#ifdef PERIODX
    !-----------------------------------
    !   Left BC
    if (coords(0).eq.0) then
       u(5,0,:,:)=u(5,nx,:,:)
    endif
    !   Right BC
    if (coords(0).eq.mpicol-1) then
       u(5,nxp1,:,:)=u(5,1,:,:)
    endif
    !-----------------------------------
#endif
#ifdef PERIODY
    !-----------------------------------
    !   bottom BC
    if (coords(1).eq.0) then
       u(5,:,0,:)= u(5,:,ny,:)
    endif
    !   top BC
    if (coords(1).eq.mpirow-1) then
       u(5,:,nyp1,:)= u(5,:,1,:)
    endif
    !-----------------------------------
#endif
#ifdef PERIODZ
    !-----------------------------------
    !   out BC
    if (coords(2).eq.0) then
       u(5,:,:,0)= u(5,:,:,nz)
    endif
    !   in BC
    if (coords(2).eq.mpirowz-1) then
       u(5,:,:,nzp1)= u(5,:,:,1)
    endif
    !-----------------------------------
#endif
    !******************************************************************
    !
#endif   !MPIP
    !
    !******************************************************************

    !******************************************************************
    !   reflecting and outflow BCs
    !-----------------------------------
    !   left
    if (coords(0).eq.0) then
       u(5,0,   0:nyp1,0:nzp1)=u(5,1 ,0:nyp1,0:nzp1)
    endif
    !-----------------------------------
    !   right
    if (coords(0).eq.mpicol-1) then
       u(5,nxp1,0:nyp1,0:nzp1)=u(5,nx,0:nyp1,0:nzp1)
    endif
    !-----------------------------------
    !   bottom
    if (coords(1).eq.0) then
       u(5,0:nxp1,0   ,0:nzp1)=u(5,0:nxp1,1 ,0:nzp1)
    endif
    !-----------------------------------
    !   top
    if (coords(1).eq.mpirow-1) then
       u(5,0:nxp1,nyp1,0:nzp1)=u(5,0:nxp1,ny,0:nzp1)
    endif
    !-----------------------------------
    !   out
    if (coords(2).eq.0) then
       u(5,0:nxp1,0:nyp1,0   )=u(5,0:nxp1,0:nyp1,1 )
    endif
    !-----------------------------------
    !   in
    if (coords(2).eq.mpirowz-1) then
       u(5,0:nxp1,0:nyp1,nzp1)=u(5,0:nxp1,0:nyp1,nz)
    endif
    !
    return
    !
  end subroutine thermal_bounds
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  Returns the length of the superstep with N inner substeps
  !--------------------------------------------------------------------
  real function superstep(N,snu)
    implicit none
    integer :: N
    real,    intent(in) :: snu  ! sqrt(nu)

    superstep=float(N)/(2.*snu) * ( (1+snu)**(2*N) - (1-snu)**(2*N) ) / &
         ( (1+snu)**(2*N) + (1-snu)**(2*N) )

    !1/( (nu-1.)*Cos(pi*(2*float(j)-1.)/(2.*float(N)) )+nu+1. )

  end function superstep
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  Returns the size of the substeps
  !--------------------------------------------------------------------
  real function substep(j,N,nu)
    implicit none
    integer, intent(in) :: j, N
    real,    intent(in) :: nu

    substep=1./( (nu-1.)*Cos(pi*float(2*j-1)/(2.*float(N)) )+nu+1. )

  end function substep
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  Returns the number of Supersteps
  !--------------------------------------------------------------------
  subroutine ST_steps(fs,Ns,fstep)
    implicit none
    real ,   intent(in) :: fs
    integer, intent(out):: Ns
    real,    intent(out):: fstep
    integer, parameter :: jmax=199
    integer :: j
    !
    do j=1,jmax
       if (superstep(j,snu) > fs) exit       
    end do

    Ns = j
    fstep = 1./superstep(Ns,snu)

  end subroutine ST_steps
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  This routine adds the heat conduction sources receives the hydro
  !  timestep in seconds, and assumes the primitives and Temp(i,j,k)
  !  arrays are updated
  !--------------------------------------------------------------------
  subroutine thermal_conduction(dt)
    implicit none
    real, intent(in) :: dt
    real :: dtcond, dts,fstep
    integer :: n,i,j,k, nsteps    
    logical :: SuperStep
    !  get the conduction timescale
    call dt_cond(dtcond)

!!$    if (dtcond < dt) then 
!!$       SuperStep=.True.
!!$       call ST_steps(dt/dtcond,Nsteps,fstep)
!!$    else
!!$       SuperStep=.False.
!!$       fstep=1.
!!$       nsteps=1
!!$    end if

    !Superstepping off
    fstep=min(dt/dtcond+1.,100.)
    nsteps=int(fstep)
    !  here i take care of the transformation from cgs to code units
    dts=dt/float(nsteps)/Psc/rsc  

    if (rank == master) print'(a,es10.3,a,i4,a)',"dt/dtcond= ", dt/dtcond, " => ", nsteps," Ssteps"

    steps : do n=1,nsteps

!!$       if (SuperStep) then
!!$          dts=dt*fstep*substep(n,Nsteps,nu)/Psc/rsc
!!$       else
!!$          dts=dt/Psc/rsc
!!$       end if

       !  show progress bar (comment to run in batch)
       !if (rank == master )call progress(n,nsteps)

       !  get the heat fluxes
       call heatfluxes()
       !  update the conserved and primitives
       do i=1,nx
          do j=1,ny
             do k=1 ,nz
                u(5,i,j,k)=u(5,i,j,k)-dts*( ( f(5,i,j,k) - f(5,i-1,j,k) )/dx &
                     +( g(5,i,j,k) - g(5,i,j-1,k) )/dy &
                     +( h(5,i,j,k) - h(5,i,j,k-1) )/dz )
             end do
          end do
       end do

       !  boundary conditions
       !  (only one layer of u(5,:,:,:) is exchanged )
       call thermal_bounds()
       !  update primitives and Temperature
       call calcprim(u, primit)

    end do steps

  end subroutine thermal_conduction
  !--------------------------------------------------------------------
end module thermal_cond
#endif
!=======================================================================
