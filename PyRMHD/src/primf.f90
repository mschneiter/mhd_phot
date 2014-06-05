!=======================================================================
!   calculates the F fluxes from the primitive variables
!   one point
!=======================================================================
subroutine primf(prim,ff)
  use parameters
  implicit none
  real,    dimension(neq), intent(in)  :: prim
  real,    dimension(neq), intent(out) :: ff
  real :: etot
  !
  etot= 0.5*prim(1)*(prim(2)**2.+prim(3)**2.+prim(4)**2.)+cv*prim(5)
  !
  ff(1) = prim(1)*prim(2)
  ff(2) = prim(1)*prim(2)*prim(2)+prim(5)
  ff(3) = prim(1)*prim(2)*prim(3)
  ff(4) = prim(1)*prim(2)*prim(4)
  ff(5) = prim(2)*(etot+prim(5))
  !
#ifdef PMHD
  ff(6)=0.0
  ff(7)=prim(2)*prim(7)-prim(6)*prim(3)
  ff(8)=prim(2)*prim(8)-prim(6)*prim(4)
#endif
#ifdef PASSIVES
  ff(neqdyn+1:neq) = prim(neqdyn+1:neq)*prim(2)
#endif
  return
end subroutine primf
!=======================================================================
