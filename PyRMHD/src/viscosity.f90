!======================================================================
!   Adds viscosity in all domain (eta is defined in parameters.f90)
!   up(n)=up(n)+eta\nabla^2(u(n))
!======================================================================
subroutine viscosity
  use parameters
  use globals
  implicit none
  integer :: i, j, k
  !
  do i=1,nx
     do j=1,ny
        do k=1,nz
           up(:,i,j,k)=up(:,i,j,k)+eta*( u(:,i+1,j,k)+u(:,i-1,j,k)       &
                                        +u(:,i,j+1,k)+u(:,i,j-1,k)       &
                                        +u(:,i,j,k+1)+u(:,i,j,k-1)       &
                                     -6.*u(:,i,j,k) )
        end do
     end do
  end do
  !
end subroutine viscosity
!======================================================================
