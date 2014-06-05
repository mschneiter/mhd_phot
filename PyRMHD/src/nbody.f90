!=======================================================================
!  N-Body module
!=======================================================================
module nbody
#ifdef NBODY
  use parameters
  implicit none
  !  N of  N-Body  
  integer :: nb
  !   position ,  velocity ans mass
  real ,allocatable :: pos(:,:), vel(:,:), mp(:)
  real :: GG   ! this is Ggrav in code units
  real :: soft= 5.*dx
  !
contains
  !--------------------------------------------------------------------
  !   Here the parameters of the module are initialized
  !   position and velociy
  !--------------------------------------------------------------------
  subroutine init_nbody()
    implicit none
    !  we need to initialize the mass position and velocity here,
    !  maybe read from file
    nb=10
    allocate( pos(3,nb) )
    allocate( vel(3,nb) )
    allocate( mp(3,nb) )

    !  fill the arrays
    ! ...

    ! change to code units
    pos(:,:)=pos(:,:)/rsc
    vel(:,:)=vel(:,:)/sqrt(vsc2)
    mp(:)=mp(:)/rhosc/(rsc**3)
    GG=Ggrav*rhosc*Rsc**2/vsc2
    !
  end subroutine init_nbody
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  computes the timestep for the particles
  subroutine nbtimestep(dt_nb)
    implicit none
    real, intent(out) :: dt_nb
    integer, i, j
    real ::rx, ry, rz, vi, vj
    
    dt_nb=1.e30
    do i=1, nb
       do j=1,i-1
          rx=pos(1,i)-pos(1,j)
          ry=pos(2,i)-pos(2,j)
          rz=pos(3,i)-pos(3,j)
          r3=sqrt(rx**2 + ry**2 + rz**2)
          vi=sqrt( vel(1,i)**2 + vel(2,i)**2 + vel(3,i)**2 + 1.0)
          vj=sqrt( vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2 + 1.0)
          !
          dt=min(dt, r3/vi)
          dt=min(dt, r3/vj)
          !
       end do
    end do
    !
  end subroutine nbtimestep
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  This surboutine calculates the acceleration on each
  !  of the  particles, at positions posp
  !--------------------------------------------------------------------
  subroutine nbaccel(posp,anb)
    implicit none
    real, intent(in) :: posp(3,nb)
    real, intent(out) :: anb(3,nb)

    integer :: i, j
    real :: rx, ry, rz, r2
    !
    !   initialize to zero
    anb(:,:)=0.
    !
    !  compute force
    do i=1,nb
       do j=1,i-1
          !
          rx=posp(1,i)-posp(1,j)
          ry=posp(2,i)-posp(2,j)
          rz=posp(3,i)-posp(3,j)
          r3=sqrt(rx**2 + ry**2 + rz**2 +soft**2)
          ! 
          !   acceleration of  i due to j
          anb(1,i) = anb(1,i) - (GG*mp(j)*rx)/(r3**3)
          anb(2,i) = anb(2,i) - (GG*mp(j)*ry)/(r3**3)
          anb(3,i) = anb(3,i) - (GG*mp(j)*rz)/(r3**3)
          !   acceletarion of j due to i
          anb(1,j) = anb(1,j) + (GG*mp(i)*rx)/(r3**3)
          anb(2,j) = anb(2,j) + (GG*mp(i)*ry)/(r3**3)
          anb(3,j) = anb(3,j) + (GG*mp(i)*rz)/(r3**3)
          !
    end do
  end do

    !
end subroutine nbaccel
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  This is a simple leap frog integrator
  !  up dates pos, vel globals to the module
  !  it takes the force on the particles an argument as well
  !  as the timestep
  !--------------------------------------------------------------------
  subroutine leapfrog(dt)
    implicit none
    real, intent(in) :: dt
    real :: posp(3,nb), aa(3,nb)
    
    do i=1, nb

       posp(:,i) = pos(:,i)+0.5*dt*vel(:,i)
       
       call accel(posp,aa)
       vel(:,i)=vel(:,i)+dt*aa(:,i)
       
       pos(:,i)=pos(:,i)+0.5*dt*vel(:,i)

    end do
    
  end subroutine leapfrog
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  This is the routine that is called from the main program
  !  to advance the position of the NBodies
  !--------------------------------------------------------------------
  begin subroutine nbody(dt_hydro)
    implicit none
    real, intent(in) :: dt_hydro
    real :: dt_nb

    call nbtimestep(dt_nb)

    if (dt_nb  < dt_hydro) then
       print*, 'NB timestep is smaller than hydro dt), stopping'
       stop
    end if

    call leapfrog(dt_hydro)



  end subroutine nbody
#endif
end module nbody
  !--------------------------------------------------------------------
