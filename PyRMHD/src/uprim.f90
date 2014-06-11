!========================================================================
!   calculates the primitive variables (rho,u,v,w,pgas) from the
!   integration variables (rho,u rho, v rho, w rho, etot)
!========================================================================
subroutine uprim(prim,uu,T)
  use parameters
  implicit none
  real,    intent(out), dimension(neq)  :: prim
  real,    intent(in),  dimension(neq)  :: uu  
  real,    intent(out)                  :: T
  real :: r, dentot
  !
  r=max(uu(1),1e-15)
  prim(1)=r
  prim(2)=uu(2)/r 
  prim(3)=uu(3)/r
  prim(4)=uu(4)/r
  !
#ifdef MHD
  !  if MHD defined the magnetic energy is included
prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2)   & 
               -0.5*  (  uu(6)**2+  uu(7)**2  +uu(8)**2) ) /cv  
#else
  ! if MHD (active) not defined
  prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2) ) /cv
#endif
  !if(prim(neqdyn).lt.0.) write(*,*) 'ay !!!'
  prim(5)=max(prim(5),1e-16)
  !
!IF PASSIVE MHD OR MHD CHOSEN
#if defined(PMHD) || defined(MHD) 
  prim(6:8) = uu(6:8)
#endif 
!
! IF EMPLOYING PASSIVE SCALARS
#ifdef PASSIVES
  prim(neqdyn+1:neq) = uu(neqdyn+1:neq)
#endif
  !
  !-----------------------------------------
#ifdef ADIABATIC
  !
  T=(prim(5)/r)*Tempsc
  !
#endif
  !-----------------------------------------
#ifdef COOLINGH
  !
  dentot=(2.*r-prim(neqdyn+1))
  dentot=max(dentot,1e-15)
  !
  T=max(1.,(prim(5)/dentot)*Tempsc)
  prim(5)=dentot*T/Tempsc
  !
#endif
  !----------------------------------------- 
#ifdef COOLINGDMC
  !
  ! assumes it is fully ionized
  r=max(r,1e-15)
  !
  T=max(1.,(prim(5)/r))*Tempsc
  prim(5)=r*T/Tempsc
  !
#endif  
  !-----------------------------------------
#ifdef COOLINGCHI
  !
  ! assumes it is fully ionized
  r=max(r,1e-15)
  !
  T=max(1.,(prim(5)/r)*Tempsc )
  prim(5)=r*T/Tempsc
  !
#endif  
  !-----------------------------------------
#ifdef COOLINGBBC
  !
  dentot= prim(neqdyn+1) + prim(neqdyn+2) + prim(neqdyn+3) + prim(neqdyn+4) &
       + prim(neqdyn+5) + prim(neqdyn+6)
  !
  T = (prim(5)/dentot/Rg)*vsc2
  !
#endif
  !-----------------------------------------
  !
return
end subroutine uprim  
!//////////////////////////////////////////////////////////////////////
