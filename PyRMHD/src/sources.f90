!=======================================================================
!  Module to commpute source terms
!=======================================================================
#ifdef SOURCE
module sources
  use parameters, only : neq, neqdyn, nxtot, nytot, nztot, &
                        Ggrav, rsc, rhosc, vsc2, nx, ny, nz
  use globals,    only : dx, dy, dz, coords
  implicit none
  
contains
  !--------------------------------------------------------------------
  !   gets the position and spherical radius calculated with respect to
  !   the center of the grid
  subroutine getpos(i,j,k,x,y,z,r)
    implicit none
    integer, intent(in)  :: i, j, k
    real,    intent(out) :: x, y, z, r
 
    x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
    y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
    z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
    
    r  = sqrt(x**2 +y**2 +z**2 )
    
  end subroutine getpos
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  Adds the gravitational force due to point particles, at this
  !  moment is fixed to two point sources (exoplanet)
  !--------------------------------------------------------------------
#ifdef GRAV
  subroutine grav_source(xc,yc,zc,pp,s)
    use star  ! this module contains the position of the planet
    implicit none
    real, intent(in)    :: xc, yc, zc
    real, intent(in)    :: pp(neq)
    real, intent(inout) :: s(neq)
    integer, parameter  :: nb=2   ! 2 particles
    real :: x(nb),y(nb),z(nb), GM(np), rad2(np)
    integer :: i

    GM(1)=0.3*Ggrav*MassS/rsc/vsc2
    GM(2)=Ggrav*MassP/rsc/vsc2

    !calculate distance from the sources
    ! star
    x(1)=xc
    y(1)=yc
    z(1)=zc
    rad2(1) = x(1)**2 +y(1)**2 + z(1)**2
    ! planet
    x(2)=xc-xp
    y(2)=yc
    z(2)=zc-zp
    rad2(2) = x(2)**2 +y(2)**2 + z(2)**2

    ! update source terms
    do i=1, nb
      ! momenta
      s(2)= s(2)-pp(1)*GM(i)*x(i)/(rad2(i)**1.5)
      s(3)= s(3)-pp(1)*GM(i)*y(i)/(rad2(i)**1.5)
      s(4)= s(4)-pp(1)*GM(i)*z(i)/(rad2(i)**1.5)
      ! energy
      s(5)= s(5)-pp(1)*GM(i)*( pp(2)*x(i) +pp(3)*y(i) +pp(4)*z(i) )  &
             /(rad2(i)**1.5 )
    end do


  end subroutine grav_source
#endif
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  Adds the radiaiton pressure force due to photo-ionization
  !--------------------------------------------------------------------
#ifdef RADPRES  
  subroutine radpress_source(i,j,k,xc,yc,zc,rc,pp,s)
#ifdef RADDIFF
  use difrad
#else
#ifdef RADDIFF_GLOBAL
  use difrad_global
#endif
#endif
    implicit none
    integer, intent(in)  :: i,j,k
    real,    intent(in)  :: xc, yc, zc, rc, pp(neq)
    real,    intent(inout) :: s(neq)
    real :: radphi
    !  the following is h/912Angstroms = h/lambda
    real :: hlambda= 7.265e-22, Frad
#ifdef RADDIFF
  radphi= ph(i,j,k) 
#else
#ifdef RADDIFF_GLOBAL
  radphi= ph(i+coords(0)*nx,j+coords(1)*ny,k+(coords(2)*nz )
#endif
#endif

  !  update the source terms
  Frad=radphi*pp(6)*hlambda*rsc/rhosc/vsc2
  !  momenta
  s(2) = s(2) + Frad*xc/rc
  s(3) = s(3) + Frad*yc/rc
  s(4) = s(4) + Frad*zc/rc
  !  energy
  s(5) = s(5)+  Frad*( xc*pp(2) + yc*pp(3) + zc*pp(4) )/rc

  end subroutine radpress_source
#endif
  !--------------------------------------------------------------------
  !  Adds terms proportional to div B in Faraday's Law,
  !  momentum equationand energy equation as propoed 
  !  in Powell et al. 1999
  !--------------------------------------------------------------------
#ifdef DIVBCORR
  subroutine divbcorr_source(i,j,k,pp,s)
   
    implicit none
    real, intent(in)  :: pp(neq)
    real, intent(out) :: s(neq)
    real              :: divB
    integer :: i, j, k

    call divergence_B(i,j,k,divB)

    ! update source terms
      ! momenta
      s(2)= s(2)-divB*pp(2) 
      s(3)= s(3)-divB*pp(3) 
      s(4)= s(4)-divB*pp(4) 

      ! energy
      s(5)= s(5)-divB*(pp(2)*pp(6)+pp(3)+pp(7)+pp(4)*pp(8))      

      ! Faraday law
      s(6)=s(6)-divB*pp(6)
      s(7)=s(7)-divB*pp(7)
      s(8)=s(8)-divB*pp(8)

  end subroutine divbcorr_source

  subroutine divergence_B(i,j,k,d)
  use globals
  implicit none
  integer, intent(in) :: i,j,k
  real, intent(out)   :: d
    
  d=  (primit(6,i,j,k)-primit(6,i-1,j,k))/dx  &
    + (primit(7,i,j,k)-primit(7,i,j-1,k))/dy  &
    + (primit(8,i,j,k)-primit(8,i,j,k-1))/dz

  end subroutine divergence_B


#endif
  !--------------------------------------------------------------------



  !--------------------------------------------------------------------
  !  Main driver, this is called from the upwind stepping (step.f90)
  !  subroutine
  !--------------------------------------------------------------------
  subroutine source(i,j,k,prim,s)
    implicit none
    integer, intent(in)  :: i, j, k
    real, intent(in)     :: prim(neq)
    real, intent(out)    :: s(neq)
    real :: x, y, z, r
    
    ! resets the source terms
    s(:) = 0.
    
    ! position with respect to the center of the grid
    call getpos( i, j, k, x, y ,z, r) 

#ifdef GRAV
    !  point source(s) gravity
    call grav_source(x,y,z,prim,s)
#endif

#ifdef RADPRES
    !  photoionization radiation pressure
    call radpress_source(i,j,k,x,y,z,r,prim,s)
#endif
#ifdef DIVBCORR
   
    !  divergence correction Powell et al. 1999
    call divbcorr_source(i,j,k,prim,s)
#endif

    
    return
  end subroutine source
  !--------------------------------------------------------------------
  
end module sources
#endif
!=======================================================================
