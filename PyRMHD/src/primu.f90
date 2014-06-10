!=======================================================================
!   calculates conserved variables from the primitives at one point
!=======================================================================
subroutine primu(prim,uu)
  use parameters
  implicit none
  real, dimension(neq), intent(in)  :: prim
  real, dimension(neq), intent(out) :: uu
  !
  uu(1) = prim(1)
  uu(2) = prim(1)*prim(2)
  uu(3) = prim(1)*prim(3)
  uu(4) = prim(1)*prim(4)
  uu(5) = 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)+cv*prim(5)
  !

#ifdef PMHD
  uu(6)=prim(6)
  uu(7)=prim(7)
  uu(8)=prim(8)
#endif

#ifdef MHD
  uu(6)=prim(6)
  uu(7)=prim(7)
  uu(8)=prim(8)
#endif

#ifdef PASSIVES
  uu(neqdyn+1:neq) = prim(neqdyn+1:neq)
#endif
!
  return
end subroutine primu
!=======================================================================
