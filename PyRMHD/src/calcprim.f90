!=======================================================================
!   calculates the primitives from the conserved vars on all domain
!=======================================================================
subroutine calcprim(u,primit)
  use parameters
#ifdef THERMAL_COND
  use thermal_cond, only : Temp
#endif
  implicit none
  real,intent(in), dimension(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :: u
  real,intent(out),dimension(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :: primit
#ifndef THERMAL_COND
  real                 :: T
#endif
  integer :: i,j,k
  !
  do i=nxmin,nxmax
     do j=nymin,nymax
        do k=nzmin,nzmax
           !
#ifdef THERMAL_COND
           call uprim(primit(:,i,j,k),u(:,i,j,k),Temp(i,j,k) )
#else
           call uprim(primit(:,i,j,k),u(:,i,j,k),T)
#endif
           !
        end do
     end do
  end do
  !
end subroutine calcprim
!=======================================================================

