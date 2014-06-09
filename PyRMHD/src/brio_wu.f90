!=======================================================================
!  Implements the BRIO-WU shock test
!=======================================================================
module BrioWu
  use parameters
  implicit none
  real  :: RSW, TSW, VSW, dsw, b0 ! star: Radius, Temperature, Velocity, Density, Magnetic Field
  real  :: RPW, TPW, VPW, dpw ! planet: Radius, Temperature, Velocity, Density
  real  :: torb, rorb ! planet: orbital period & distance
  real  :: MassS, MassP
  real  :: xp, yp, zp  ! position of the planet (respect to the star)
contains

  subroutine init_bw()
    implicit none
    real :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, pL, pR
    real :: BxL, BxR, ByL, ByR, BzL, BzR

    rhoL = 1.
    rhoR = 0.125
    vxL  = 0.
    vxR  = 0.
    vyL  = 0.
    vyR  = 0.
    vzL  = 0.
    vzR  = 0.
    pL   = 1.
    pR   = .1
    BxL  = 0.
    BxR  = 0.
    ByL  = 1.
    ByR  = -1.
    BzL  = 0.
    BzR  = 0.

  end subroutine init_bw

  !--------------------------------------------------------------------

  subroutine impose_bw(u,time)
    use globals, only : coords, dx, dy, dz, rank
    use parameters, only : nxtot
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: time
    real :: x !, y, z, xpl, ypl, zpl
! BRIO-WU 
    real :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, pL, pR
    real :: BxL, BxR, ByL, ByR, BzL, BzR

    integer ::  i,j,k

    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax

             ! Position measured from the centre of the grid (star)
             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx

             ! IF LEFFT STATE
             if( x <= nxtot/2) then

                !   total density and momena
                u(1,i,j,k) = rhoL
                u(2,i,j,k) = rhoL*vxL
                u(3,i,j,k) = rhoL*vyL
                u(4,i,j,k) = rhoL*vzL
                !   energy
                u(5,i,j,k)=0.5*rhoL*(vxL**2+vyL**2+vzL**2)+cv*pL
#ifdef MHD
                !  add magnetic energy
                u(5,i,j,k) =  u(5,i,j,k)+0.5*(BxL**2+ByL**2+BzL**2)

                u(6,i,j,k) =  BxL
                u(7,i,j,k) =  ByL
                u(8,i,j,k) =  BZL
#endif
             else !RIGHT STATE
                !   total density and momena
                u(1,i,j,k) = rhoR
                u(2,i,j,k) = rhoR*vxR
                u(3,i,j,k) = rhoR*vyR
                u(4,i,j,k) = rhoR*vzR
                !   energy
                u(5,i,j,k)=0.5*rhoR*(vxR**2+vyR**2+vzR**2)+cv*pR

#ifdef MHD
                !  add magnetic energy
                u(5,i,j,k) =  u(5,i,j,k)+0.5*(BxR**2+ByR**2+BzR**2)

                u(6,i,j,k) =  BxR
                u(7,i,j,k) =  ByR
                u(8,i,j,k) =  BZR

#endif
             end if
                
          end do
       end do
    end do

  end subroutine impose_bw
  !--------------------------------------------------------------------

end module  BrioWu
