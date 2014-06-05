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
  prim(5)=( ( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2) ) /cv )
  !if(prim(neqdin).lt.0.) write(*,*) 'ay !!!'
  prim(5)=max(prim(5),1e-16)
  !
#ifdef PMHD
  prim(6) = uu(6)
  prim(7) = uu(7)
  prim(8) = uu(8)            
#endif 


#ifdef MHD
  prim(6) = uu(6)
  prim(7) = uu(7)
  prim(8) = uu(8)            
#endif 

!
#ifdef PASSIVES
  prim(neqdin+1:neq) = uu(neqdin+1:neq)
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
  dentot=(2.*r-prim(neqdin+1))
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
  dentot= prim(neqdin+1) + prim(neqdin+2) + prim(neqdin+3) + prim(neqdin+4) &
       + prim(neqdin+5) + prim(neqdin+6)
  !
  T = (prim(5)/dentot/Rg)*vsc2
  !
#endif
  !-----------------------------------------
  !
return
end subroutine uprim  
!//////////////////////////////////////////////////////////////////////
