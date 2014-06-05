!=======================================================================
!   integration from t to t+dt with an approximate Riemann solver
!=======================================================================
subroutine tstep(time,dt)
  use parameters
  use globals
#ifdef RADDIFF
  use difrad
#endif
#ifdef RADDIFF_GLOBAL
  use difrad_global
#endif
#ifdef THERMAL_COND
  use thermal_cond
#endif
  implicit none
  real, intent(in):: dt, time
  real            :: dtm
  !
  !---------------------------------------------------------------------
  !  1st half timestep
  !---------------------------------------------------------------------
  !
  dtm=dt/2.
  !
  !   calculate the fluxes using the primitives
  !   (piecewise constant)
  !   ---> f_i^n=Riemann[p_i^n,p_{i+1}^n]
#ifdef HLL
  call hllfluxes(1)
#endif
#ifdef HLLC
  call hllcfluxes(1)
#endif
  !   upwind timestep
  !   up=u^{n+1/2}=u_i^n-(dt/2dx)[f_i^n - f_{i-1}^n]
  call step(dtm)
  !
  !   add viscosity
  !call viscosity
  !
  !---------------------------------------------------------------------
  !  2nd half timestep
  !---------------------------------------------------------------------
  !
  !   boundaries in up and  primitives up ---> primit
  call boundaryII(time,dt)
  call calcprim(up,primit)
  !
  !   calculate the fluxes using the primitives
  !   with linear reconstruction (piecewise linear)
  !   p_{R/L}=p_{i(+/-)1/2}^{n+1/2}
  !   ---> f_{i+1/2}^(n+1/2)=Riemann[p_{R,i}^{n+1/2},p_{L+1}^{n+1/2}]
#ifdef HLL
  call hllfluxes(2)
#endif
#ifdef HLLC
  call hllcfluxes(2)
#endif
  !   up=u^{n+1}=u_i^n-(dt/dx)[f_{i+1/2}^{n+1/2}-f_{i-1/2}^{n+1/2}]
   call step(dt)
  !
  !   add viscosity
  call viscosity
  !
  !--------------------------------------------------------------------
  !   copy the up's on the u's
  u=up
  !
#ifdef RADDIFF
  call diffuse_rad()
#endif
#ifdef RADDIFF_GLOBAL
  call diffuse_rad()
#endif
  !****************************************************
  !   apply cooling/heating
  !****************************************************
#ifdef COOLINGH
  !   add cooling to the conserved variables
  call coolingh(dt*tsc)
#endif
  !****************************************************
#ifdef COOLINGDMC
  !   the primitives are updated in the cooling routine
  call coolingdmc(dt*tsc)
#endif
#ifdef COOLINGCHI
  !   the primitives are updated in the cooling routine
  call coolingchi(dt*tsc)
#endif
 !****************************************************
  !   BBC cooling not implemented yet
  !****************************************************
  !
  !----------------------------------------------------
  !   boundary contiditions on u
  call boundaryI(time,dt)
  !----------------------------------------------------
  !
  !  update the primitives with u
  call calcprim(u, primit)
  !
  !----------------------------------------------------
  !  Thermal conduction
#ifdef THERMAL_COND
  call thermal_conduction(dt*tsc)
#endif
  !----------------------------------------------------
  !
end subroutine tstep
!=======================================================================
