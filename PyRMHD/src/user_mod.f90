!=======================================================================
!  User imput module
!  This is an attempt to have all input neede from user in a single file
!  This module should load additional modules (i.e. star, jet, sn), to
!  impose initial and boundary conditions (such as sources)
!=======================================================================
module user_mod
  use parameters
  !  here we should load additional modules
  !use star
  !use jet
  use BrioWu
  implicit none
  !  Variables global to the module can be here
  ! real :: variablename

contains

  !--------------------------------------------------------------------
  !   Here the parameters of the User Defined module are initialized, 
  !   and scaled to code units
  !   This routine is alway called initmain.f90
  !   It has to be present, even if empty
  !--------------------------------------------------------------------
  subroutine init_user_mod()
    implicit none
      
    !  initialize modules loaded by user
    call init_bw()

  end subroutine init_user_mod
  !--------------------------------------------------------------------



  !--------------------------------------------------------------------
  !   Here the domain is initialized at t=0
  !   In general takes U and time as argument, if time is not needed,
  !   just don't use it...
  !--------------------------------------------------------------------
  subroutine initial_conditions(u,time)
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in) :: time

    
    !   Fills the domain w/the setup of the Brio-Wu shock
    call impose_bw(u)
    
  end subroutine initial_conditions
  !--------------------------------------------------------------------
 


  !--------------------------------------------------------------------
  !   User Defined Boundary conditions 
  !--------------------------------------------------------------------
#ifdef OTHERB
  subroutine impose_user_bc(u,time)
   use globals, only : coords, dx, dy, dz, rank
   implicit none
   real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
   real, intent(in) :: time
   !integer ::  i,j,k
   
   !  We can impose here the exoplanet, has to load star module at the
   !  top
   !call impose_star(u,time)
   
  end subroutine impose_user_bc
#endif
  !--------------------------------------------------------------------

end module user_mod
