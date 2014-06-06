!=======================================================================
!   calculates the sound speed or fast magnetosonic speed
!=======================================================================
module sound
use parameters
implicit none

contains

!=======================================================================
!   calculates the sound speed for HD and PMHD
!=======================================================================
subroutine csound(p,d,cs)
!  use parameters
  implicit none
  real, intent(in) :: p, d
  real, intent(out) ::cs
  !
  cs=sqrt(gamma*p/d)
  !
 ! return
end subroutine csound
!=======================================================================

!=======================================================================
!   calculates the magnetosonic speed for MHD
!=======================================================================
subroutine cfast(p,d,bx,by,bz,cfx,cfy,cfz)
!  use parameters
  implicit none
  real, intent(in) :: p, d, bx, by, bz
  real, intent(out) ::cfx,cfy,cfz
  real :: b2
  !
  b2=bx*bx+by*by+bz*bz
  cfx=sqrt(0.5*((gamma*p+b2)/d+sqrt((gamma*p+b2)/d)**2-4.*gamma*p*bx*bx/d**2))
  cfy=sqrt(0.5*((gamma*p+b2)/d+sqrt((gamma*p+b2)/d)**2-4.*gamma*p*by*by/d**2))
  cfz=sqrt(0.5*((gamma*p+b2)/d+sqrt((gamma*p+b2)/d)**2-4.*gamma*p*bz*bz/d**2))                
  !
  !  return
end subroutine cfast

end module sound
!=======================================================================
