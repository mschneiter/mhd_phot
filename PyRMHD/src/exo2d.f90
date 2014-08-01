!=======================================================================
!  Implements the ORZAG-TANG VORTEX test
!=======================================================================
module Exo2D
  use parameters
  implicit none
!  real  :: RSW, TSW, VSW, dsw, b0 ! star: Radius, Temperature, Velocity, Density, Magnetic Field
  real  :: RPW, TPW, VPW, dpw ! planet: Radius, Temperature, Velocity, Density
  real  :: torb, rorb ! planet: orbital period & distance
  real  :: MassS, MassP
  real  :: xp, yp, zp  ! position of the planet (respect to the star)


contains

  subroutine init_exo2d()
    implicit none
    
!PLANETARY PARAMETERS
    real :: ampdot ! mdot_planet (g/s)

    !PLANETS PARAMETERS              (Pascal Tremblin & Eugene Chiang 2013)
    MassP=0.67*mjup
    AMPDOT=1.6E11      !3.E10          ! Planetary Mass Loss rate (g/s) !Linksys. et al 2010
    TPW  = 7.E3        !1.E4           ! Planets temperature
    !RPW  = 0.002002*AU!50.E12         ! Planetary wind radius (cm)   
    RPW  =  4.*Rexo                    ! Planetary wind radius (cm)
    vpw  = 12.e5    !10.E5             ! Planets wind velocity (cm/s)
    dpw=((AMDOT/RPW)/(4*pi*RPW*VPW))   ! Planetary density
 
    ! change to code units
    dpw=dpw/rhosc
    vpw=vpw/sqrt(vsc2)
    Tpw=Tpw/Tempsc
    Rpw=Rpw/rsc

    !ORBITAL PARAMETERS
    rorb=.047*AU!0.47**AU
    torb=3.52*day

    ! change to code units
    rorb=rorb/rsc 
    torb=torb/tsc

    !  initial position
    xp=Rorb*cos(-25.*pi/180.)
    yp=0.
    zp=Rorb*sin(-25.*pi/180.)
        
  end subroutine init_exo2d

  !--------------------------------------------------------------------
  ! Initial conditions for the BRIO-WU shock tube
  subroutine impose_exo2d(u,time)
    use globals, only : coords, dx, dy, rank
    use parameters, only : nxtot
    implicit none

!!    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
!!    real :: x, y
!!    integer ::  i,j,k
!!
    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax

             ! Position measured from the bottom corner of the grid
             x=(float(i+coords(0)*nx)+0.5)*dx*rsc !WHY +0.5?
             y=(float(j+coords(1)*ny)+0.5)*dy*rsc !WHY +0.5?
                                              !ANS: Position of cell center
             

             !IMPOSE STELLAR WIND FROM THE -X-SIDE
             if ( x <= 10.) then
                vx = 130D5
                vy = 0d0
                vz = 0d0
                bx = 0d0
                by = 10d0/sqrt(4*pi)
                bz = 0.
                              

                !   total density and momena
                u(1,i,j,k) = rho
                u(2,i,j,k) = rho*vx
                u(3,i,j,k) = rho*vy
                u(4,i,j,k) = rho*vz
                !   total energy
                u(5,i,j,k)=0.5*rho*(vx**2+vy**2+vz**2)+cv*p  &
                               +0.5*(Bx**2+By**2+Bz**2)
                u(6,i,j,k) =  Bx
                u(7,i,j,k) =  By
                u(8,i,j,k) =  BZ

          end do
       end do
    end do

  end subroutine impose_ot
  !--------------------------------------------------------------------

end module OrzagTang
