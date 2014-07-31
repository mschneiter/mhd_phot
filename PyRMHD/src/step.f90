!=======================================================================
!   Performs the upwind time step 
!=======================================================================
subroutine step(dt)
  use parameters
  use globals
#ifdef SOURCE
  use sources
#endif
  implicit none
#ifdef SOURCE
  real :: s(neq)
#endif
  real, intent(in) :: dt
  integer :: i, j, k
  real :: dtdx, dtdy, dtdz
  !
  dtdx=dt/dx
  dtdy=dt/dy
  dtdz=dt/dz
  !
  !------------------------------------------------------------------  
  do i=1,nx
     do j=1,ny
        do k=1,nz
           !
           up(:,i,j,k)=u(:,i,j,k)-dtdx*(f(:,i,j,k)-f(:,i-1,j,k))    &
                                 -dtdy*(g(:,i,j,k)-g(:,i,j-1,k))    &
                                 -dtdz*(h(:,i,j,k)-h(:,i,j,k-1))

#ifdef SOURCE
           call source(i,j,k,primit(:,i,j,k),s)
           up(:,i,j,k)= up(:,i,j,k)+dt*s(:)
#endif
           !
        end do
     end do
  end do
  !
end subroutine step


!=======================================================================


! Calculates the divergence of the magnetic field

!! subroutine calcdivb(i,j,k,primit(:,i,j,k),divb(i,j,k))
!!   use globals, only : dx, dy, dz
!!   implicit none
!!   real, dimension(nx,ny,nz) :: divb
!!  
!!  divb(i,j,k) =    (primit(6,i,j,k)-primit(6,i-1,j,k))/dx  &
!!                 + (primit(7,i,j,k)-primit(7,i,j-1,k))/dy  &
!!                 + (primit(8,i,j,k)-primit(8,i,j-1,k))/dz
!!   
!!  
!!  end subroutine divb

 
