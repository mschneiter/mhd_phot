!=======================================================================
!   swapping routine (in y-direction)
!=======================================================================
subroutine swapy(prim,neq)
  implicit none
  real, intent(inout), dimension(neq) :: prim
  integer, intent(in) :: neq
  real :: aux
  !
  aux=prim(2)
  prim(2)=prim(3)
  prim(3)=aux
  !
#ifdef PMHD
  aux=prim(6)
  prim(6)=prim(7)
  prim(7)=aux
#endif 

  return
end subroutine swapy
!=======================================================================

!=======================================================================
!   swapping routine (in z-direction)
!=======================================================================
subroutine swapz(prim,neq)
  implicit none
  real, intent(inout), dimension(neq) :: prim
  integer, intent(in) :: neq
  real :: aux
  !
  aux=prim(2)
  prim(2)=prim(4)
  prim(4)=aux
  !
#ifdef PMHD
  aux=prim(6)
  prim(6)=prim(8)
  prim(8)=aux
#endif

  return
end subroutine swapz
!=======================================================================
