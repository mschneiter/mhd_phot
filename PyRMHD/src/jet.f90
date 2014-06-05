!=======================================================================
!  Jet source module
!=======================================================================
module jet
  use parameters
  implicit none
  real, save :: Rj, Lj, denj, Tempj, vj0, dvj, tau, omega 
  real, save :: posj(3)
  !  the direction can be obtained with the following parameters
  !  alpha is the angle with respect to z at t=0
  !  the angle can be adjusted by rotating it by a precesion
  !  angle (omegaP x t)
  real, save :: alpha, omegaP

contains
  !--------------------------------------------------------------------
  !   Here the parameters of the jet are initialized, and scaled to
  !   code units
  !--------------------------------------------------------------------
  subroutine init_jet()
    implicit none
    
    Rj    = 3.e16/rsc   !  jet radius
    Lj    = 3.e16/rsc   !  jet length
    
    !  jet position
    posj(1)= 5.e17/rsc
    posj(2)= 5.e17/rsc
    posj(3)= 0.e17/rsc
    
    !  jet orientation parameters
    alpha =5.*pi/180.
    omegaP=2.*pi/(400.*yr/tsc)

    denj  = 100.                    !  density
    Tempj = 1000.                  !  jet temperature
    vj0   = 500.e5/sqrt(vsc2)      !  mean velocity
    dVj   = 250.e5/sqrt(vsc2)      !  amplitude of variability
    tau   = 150.*yr/tsc            !  period of variability
    omega = 2.*pi/tau              !
    

  end subroutine init_jet
  !--------------------------------------------------------------------
  subroutine impose_jet(u,time)
    use globals, only : coords, dx, dy, dz, rank
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: time
    real :: omegat, x, y, z, rad, xp, yp, zp, rx, ry, rz, vjet
    !  precesion opening angle (or initial direction, repect to the z axis)
    real :: sina, cosa
    real :: coso, sino
    integer ::  i,j,k


    sina= sin(alpha)
    cosa= cos(alpha)
    
    sino= sin(omegaP*time)
    coso= cos(omegaP*time)

    omegat = omega*time   ! for jet variability

    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax
           
             !   measured from the corner of the computational mesh
             x=(float(i+coords(0)*nx)+0.5)*dx
             y=(float(j+coords(1)*ny)+0.5)*dy
             z=(float(k+coords(2)*nz)+0.5)*dz

             xp=x-posj(1)
             yp=y-posj(2)
             zp=z-posj(3)

             rx= xp*coso    -yp*sino
             ry= xp*cosa*sino+yp*cosa*coso -zp*sina
             rz= xp*sina*sino+yp*sina*coso +zp*cosa
             
             rad=sqrt(rx**2+ry**2)          

             !if( (j.eq.0).and.(i.eq.0).and.(rank.eq.0)) print*,k,z,zp

             if( (abs(rz).le.Lj).and.(rad.le.(Rj)) ) then

                !  inside the jet source
                vjet= vj0+dvj*sin(omegat)
                vjet=sign(vjet,rz)
                !
                !   total density and momenta
                u(1,i,j,k) = denj
                u(2,i,j,k) = denj*vjet*sina*coso
                u(3,i,j,k) = denj*vjet*sina*sino
                u(4,i,j,k) = denj*vjet*cosa
                !   energy
                u(5,i,j,k)=0.5*denj*vjet**2+cv*denj*1.0001*Tempj/Tempsc
                !  density of neutrals
                u(6,i,j,k)= 0.9999*denj
                !   passive scalar
                u(7,i,j,k)= denj
             endif

!!$    real :: rstar=1e16 !!!!!!!!!!!!!!!!!!!!
!!$             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
!!$             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
!!$             z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
!!$             
!!$             rad=sqrt(x**2+y**2+z**2)*rsc
!!$             
!!$             if (rad.le.rstar) then
!!$                !   total density and momenta
!!$                u(1,i,j,k)= nenv
!!$                u(2,i,j,k)= 0.
!!$                u(3,i,j,k)= 0.
!!$                u(4,i,j,k)= 0.
!!$                !   energy
!!$                u(5,i,j,k)=cv*1.9999*nenv*1000./Tempsc
!!$                !  density of neutrals
!!$                u(6,i,j,k)= 0.0001*nenv
!!$                !   passive scalar
!!$                u(7,i,j,k)= -nenv
!!$             end if
             
          end do
       end do
    end do
    
    !call mpi_finalize(i)
    !stop

  end subroutine impose_jet
  !--------------------------------------------------------------------
end module jet
