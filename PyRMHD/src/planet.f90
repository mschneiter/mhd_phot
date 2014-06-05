!=======================================================================
!  Planet source module (moving wind)
!=======================================================================
module planet
  use parameters
!    use parameters, only : yr, msun, rsc, T0, vsc2, rhosc 
  implicit none
!  real  :: RSW, TSW, VSW, dsw ! star Radius, Temperature & Velocity
  real  :: RPW, TPW, VPW, DPW ! star Radius, Temperature & Velocity

contains
  !--------------------------------------------------------------------
  !   Here the parameters of the Star are initialized, and scaled to
  !   code units
  !--------------------------------------------------------------------
  subroutine init_planet()
    implicit none
    real :: amdot,m_jup, rp, rpau ! mdot_planet (MSUN/yr)
    real :: rorb
!    r_jup=71500E3!Jupiter's mass (cm)
!    rp=1.4*r_jup
    AMDOT=2.E10                       ! Planetary Mass Loss rate (g/s)
    RPW  = 1.E12                      ! Planetary wind radius (cm)   
    TPW  = 1.E4                       ! Planets temperature
    vpw  = 10.e5                      ! Planets wind velocity (cm/vsc)
    dpw=((AMDOT/RSW)/(4*pi*RSW*VSW))  ! Planetary density
    rorb=0.047*AU                     ! Orbital radius (cm)
    ! change to code units
    dpw=dpw/rhosc
    vpw=vpw/sqrt(vsc2)
    Tpw=Tpw/Tempsc
    Rpw=Rpw/rsc
    rorb=rorb/rsc
    !print*,dsw,vsw
  end subroutine init_planet
  !--------------------------------------------------------------------
  subroutine impose_planet(u,t)
    use globals, only : coords, dx, dy, dz, rank
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: t
    real :: x, y, z
    real :: velx, vely, velz, rads, dens, rorb
    integer ::  i,j,k
    
    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax
           
             ! Position measured from the centre of the planet
             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
             z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

             ! Distance from the centre of the planet
             rads=sqrt(x**2+y**2+z**2)          
             if(rads.eq.0) rads=dx*0.10
             if( rads.le.rsw) then

               VelX=VPW*X/RADS
               VelY=VPW*Y/RADS
               VelZ=VPW*Z/RADS
               DENS=DSW!*RSW**2/RADS**2
                !   total density and momenta
                u(1,i,j,k) = dens
                u(2,i,j,k) = dens*velx
                u(3,i,j,k) = dens*vely
                u(4,i,j,k) = dens*velz
                !   energy
                u(5,i,j,k)=0.5*dens*vsw**2 &
                     + cv*dens*1.9999*Tsw
                !  density of neutrals
                u(neqdyn+1,i,j,k)= 0.0001*dens
                !   passive scalar
                u(neqdyn+2,i,j,k)= dens
             endif             
          end do
       end do
    end do

    !call mpi_finalize(i)
    !stop

  end subroutine impose_star
  !--------------------------------------------------------------------
end module star
