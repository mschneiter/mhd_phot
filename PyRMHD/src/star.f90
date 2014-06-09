!=======================================================================
!  Star source module
!=======================================================================
module star
  use parameters
  implicit none
  real  :: RSW, TSW, VSW, dsw, b0 ! star: Radius, Temperature, Velocity, Density, Magnetic Field
  real  :: RPW, TPW, VPW, dpw ! planet: Radius, Temperature, Velocity, Density
  real  :: torb, rorb ! planet: orbital period & distance
  real  :: MassS, MassP
  real  :: xp, yp, zp  ! position of the planet (respect to the star)
contains
!!$  !--------------------------------------------------------------------
!!$  !   Here the parameters of the Star are initialized, and scaled to
!!$  !   code units
!!$  !--------------------------------------------------------------------
!!$  subroutine init_star()
!!$    implicit none
!!$    real :: amdot  ! mdot_star (MSUN/yr)
!!$    real :: ampdot ! mdot_planet (g/s)
!!$    real :: bsw, bsc
!!$
!!$    !STARS PARAMETERS
!!$    MassS = 1.1*msun
!!$    AMDOT = 9.E-14*msun/yr              ! Stellar Mass Loss rate (g s^-1)
!!$    TSW   = 3.E6                        ! Stellar temperature (K)
!!$    !  Stellar wind, imposed at the 1.5x  sonic point (cm)
!!$    RSW   = 1.5*Ggrav*MassS/2./(Rg*Tsw/0.6)
!!$    vsw   = 372.e5                        ! Stellar wind velocity (cm/s)
!!$    dsw   =((AMDOT/RSW)/(4*pi*RSW*VSW))   ! Stellar density @RS (g cm^-3)
!!$    bsw    =1.0e-4                         ! Stellar magnetic field (g)
!!$
!!$    ! change to code units
!!$    dsw=dsw/rhosc
!!$    vsw=vsw/sqrt(vsc2)
!!$    Tsw=Tsw/Tempsc
!!$    Rsw=Rsw/rsc
!!$    bsw=bsw/bsc 
!!$
!!$    !PLANETS PARAMETERS
!!$    MassP=0.67*mjup
!!$    AMPDOT=3.E10                       ! Planetary Mass Loss rate (g/s)
!!$    TPW  = 1.E4                        ! Planets temperature
!!$    !RPW  = 0.002002*AU!50.E12         ! Planetary wind radius (cm)   
!!$    RPW  = 3.*9.e9  *4.                ! Planetary wind radius (cm)
!!$    vpw  = 10.e5                       ! Planets wind velocity (cm/s)
!!$    dpw=((AMDOT/RPW)/(4*pi*RPW*VPW))   ! Planetary density
!!$ 
!!$    ! change to code units
!!$    dpw=dpw/rhosc
!!$    vpw=vpw/sqrt(vsc2)
!!$    Tpw=Tpw/Tempsc
!!$    Rpw=Rpw/rsc
!!$
!!$    !ORBITAL PARAMETERS
!!$    rorb=.047*AU!0.47**AU
!!$    torb=3.52*day
!!$
!!$    ! change to code units
!!$    rorb=rorb/rsc 
!!$    torb=torb/tsc
!!$
!!$    !  initial position
!!$    xp=Rorb*cos(-25.*pi/180.)
!!$    yp=0.
!!$    zp=Rorb*sin(-25.*pi/180.)
!!$
!!$  end subroutine init_star
!!$
!!$  !--------------------------------------------------------------------
!!$
!!$  subroutine impose_star(u,time)
!!$    use globals, only : coords, dx, dy, dz, rank
!!$    implicit none
!!$    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
!!$    real, intent (in)   :: time
!!$    real :: x, y, z, xpl, ypl, zpl
!!$    real :: velx, vely, velz, rads, dens, radp, phi
!!$    real :: vxorb, vyorb, vzorb, omega, cpi
!!$    integer ::  i,j,k
!!$
!!$    phi=-25.*pi/180.
!!$
!!$    omega=2.*pi/TORB
!!$
!!$    XP=Rorb*COS(omega*TIME+phi)
!!$    ZP=Rorb*SIN(omega*TIME+phi)
!!$
!!$    vxorb=-omega*Rorb*sin(omega*TIME+phi)
!!$    vzorb= omega*Rorb*cos(omega*TIME+phi)
!!$    vyorb=0.
!!$
!!$    !IF(RANK == 0) print*,'orbit ',xp,zp,vxorb*sqrt(vsc2)/1E5,vzorb*sqrt(vsc2)/1E5, omega*Rorb*sqrt(vsc2)/1E5
!!$
!!$    do i=nxmin,nxmax
!!$       do j=nymin,nymax
!!$          do k=nzmin,nzmax
!!$
!!$             ! Position measured from the centre of the grid (star)
!!$             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
!!$             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
!!$             z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
!!$             ! Position measured from the centre of the planet
!!$             xpl=x-xp
!!$             ypl=y
!!$             zpl=z-zp
!!$
!!$             ! Distance from the centre of the star
!!$             rads=sqrt(x**2+y**2+z**2)
!!$             cpi=b0*(rsw/(rads+1.e-30))**3/(2.*(rads+1.e-30)**2)
!!$ 
!!$             ! Distance from the centre of the planet
!!$             radp=sqrt(xpl**2+ypl**2+zpl**2)         
!!$             
!!$             ! IF INSIDE THE STAR
!!$             if( rads <= rsw) then
!!$                if(rads == 0.) rads=dx*0.10
!!$                
!!$                VelX=VSW*X/RADS
!!$                VelY=VSW*Y/RADS
!!$                VelZ=VSW*Z/RADS
!!$                DENS=DSW*RSW**2/RADS**2
!!$                !   total density and momena
!!$                u(1,i,j,k) = dens
!!$                u(2,i,j,k) = dens*velx
!!$                u(3,i,j,k) = dens*vely
!!$                u(4,i,j,k) = dens*velz
!!$                !   energy
!!$                u(5,i,j,k)=0.5*dens*vsw**2 &
!!$                     + cv*dens*2.*Tsw
!!$#ifdef PMHD
!!$                u(6,i,j,k) =  3.*y*x*cpi
!!$                u(7,i,j,k) = (3.*y**2-rads**2)*cpi
!!$                u(8,i,j,k) =  3.*y*z*cpi
!!$#endif
!!$
!!$                !  density of neutrals
!!$                u(neqdin+1,i,j,k)= 0.*dens
!!$                !   passive scalar
!!$                u(neqdin+2,i,j,k)= dens
!!$                
!!$                ! IF INSIDE THE PLANET
!!$             else if(radp <= rpw) then
!!$                if(radp == 0.) radp=dx*0.10
!!$                
!!$                VelX=VXORB+VPW*XPL/RADP
!!$                VelY=VYORB+VPW*YPL/RADP
!!$                VelZ=VZORB+VPW*ZPL/RADP
!!$                DENS=DPW*RPW**2/RADP**2
!!$                !   total density and momenta
!!$                u(1,i,j,k) = dens
!!$                u(2,i,j,k) = dens*velx
!!$                u(3,i,j,k) = dens*vely
!!$                u(4,i,j,k) = dens*velz
!!$                !   energy
!!$                u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
!!$                     + cv*dens* 1.0001*Tpw
!!$#ifdef PMHD
!!$                u(6,i,j,k) = 0.0
!!$                u(7,i,j,k) = 0.0
!!$                u(8,i,j,k) = 0.0
!!$#endif
!!$                !  density of neutrals
!!$                u(neqdin+1,i,j,k)=0.9999*dens
!!$                !   passive scalar
!!$                u(neqdin+2,i,j,k)= -dens
!!$                
!!$             end if
!!$                
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    !call mpi_finalize(i)
!!$    !stop
!!$
!!$  end subroutine impose_star
!!$  !--------------------------------------------------------------------

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%% BRIO-WU TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine init_star()
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

  end subroutine init_star

  !--------------------------------------------------------------------

  subroutine impose_star(u,time)
    use globals, only : coords, dx, dy, dz, rank
    use parameters, only : nxtot
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: time
    real :: x !, y, z, xpl, ypl, zpl
! BRIO-WU 
    real :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, pL, pR
    real :: BxL, BxR, ByL, ByR, BzL, BzR

!    real :: velx, vely, velz, rads, dens, radp, phi
!    real :: vxorb, vyorb, vzorb, omega, cpi
    integer ::  i,j,k
!    integer :: disc

    !IF(RANK == 0) print*,'orbit ',xp,zp,vxorb*sqrt(vsc2)/1E5,vzorb*sqrt(vsc2)/1E5, omega*Rorb*sqrt(vsc2)/1E5
    
!    disc = nxtot/2


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
#ifdef PMHD
                u(6,i,j,k) =  0.
                u(7,i,j,k) =  1.
                u(8,i,j,k) =  0.
#endif
             else !RIGHT STATE
                !   total density and momena
                u(1,i,j,k) = rhoR
                u(2,i,j,k) = rhoR*vxR
                u(3,i,j,k) = rhoR*vyR
                u(4,i,j,k) = rhoR*vzR
                !   energy
                u(5,i,j,k)=0.5*rhoR*(vxR**2+vyR**2+vzR**2)+cv*pR

#ifdef PMHD
                u(6,i,j,k) =  0.
                u(7,i,j,k) =  1.
                u(8,i,j,k) =  0.
#endif
             end if
                
          end do
       end do
    end do

    !call mpi_finalize(i)
    !stop

  end subroutine impose_star
  !--------------------------------------------------------------------








end module star
