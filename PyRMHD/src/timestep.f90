!=======================================================================
!   calculates the timestep as permited by the CFL criterion
!   Courant number set in parameters.f90
!=======================================================================
subroutine timestep(dt)
  use parameters
  use globals
  use sound
#ifdef OTHERB
  use star
#endif
  implicit none
  real, intent(out) ::dt
  real              :: dtp  
#ifdef MHD
  real              :: cx, cy, cz
#else
  real              :: c
#endif
  integer :: i, j, k, err
  !
  dtp=1.e30
  do i=1,nx
     do j=1,ny
        do k=1,nz
           !
#ifdef MHD
           call cfast(primit(5,i,j,k),primit(1,i,j,k),&
                primit(6,i,j,k), primit(7,i,j,k), primit(8,i,j,k), &
                cx,cy,cz)
           dtp=min(dtp,dx/(abs(primit(2,i,j,k))+cx))  
           dtp=min(dtp,dy/(abs(primit(3,i,j,k))+cy))
           dtp=min(dtp,dz/(abs(primit(4,i,j,k))+cz))
#else
           call csound(primit(5,i,j,k),primit(1,i,j,k),c)

           dtp=min(dtp,dx/(abs(primit(2,i,j,k))+c))  
           dtp=min(dtp,dy/(abs(primit(3,i,j,k))+c))
           dtp=min(dtp,dz/(abs(primit(4,i,j,k))+c))
#endif
!!   WE HAVE TO USE 
!!           dtp=min(dtp,dx/(abs(primit(2,i,j,k))+cs))  
!!        !
        end do
     end do
  end do
  dtp=cfl*dtp
#ifdef MPIP
  call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
  dt=dtp
#endif
  !
  return
end subroutine timestep
!=======================================================================
