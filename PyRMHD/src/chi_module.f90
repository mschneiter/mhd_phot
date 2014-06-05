!========================================================================
!   CHIANTI equilibrium cooling module
!   The table is filled at initmain.f90
!   cooldchi(T) interpolates the cooling table and returns the coefficient
!========================================================================
module chi_module
#ifdef COOLINGCHI
  !use globals
  implicit none
  real (kind=8), dimension(2,41) :: cooltab
  !
contains
  function coolchi(T)
    !use globals
    !implicit none
    real , intent(in) :: T
    integer           :: if1
    real (kind=8)     :: coolchi, T0, T1, C0, C1
    !
    if(T.gt.1e8) then
       coolchi=0.21D-26*Sqrt(dble(T))
    else
       if1=int(log10(T)*10)-39
       T0=cooltab(1,if1)
       c0=cooltab(2,if1)
       T1=cooltab(1,if1+1)
       c1=cooltab(2,if1+1)
       coolchi=(c1-c0)*(dble(T)-T0)/(T1-T0)+c0
    end if
    !
  end function coolchi
  !
#endif
end module chi_module
!========================================================================
