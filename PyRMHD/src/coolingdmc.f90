!======================================================================
!   Cooling routine foro coronal equilibrium (Dalgarno & Mc Cray 1972)
!   The cooling is turned of bellow 10^4 K
!   And the primitives are updated in this routine as well
!======================================================================
subroutine coolingdmc(dt)
#ifdef COOLINGDMC
  use parameters
  use globals
  use dmc_module
  implicit none
  real,    intent(in)  :: dt
  real                 :: T ,Eth0, dens
  real, parameter :: Tmin=10000.
  real (kind=8)        :: ALOSS, Ce
  integer :: i, j, k
  !
  do i=1,nx
     do j=1,ny
        do k=1,nz
           !
           !   get the primitives (and T)
           call uprim(primit(:,i,j,k),u(:,i,j,k),T)
           !
           if(T.gt.Tmin) then
              !
              Eth0=cv*primit(neqdyn,i,j,k)
              !
              Aloss=cooldmc(T)
              dens=primit(1,i,j,k)
              Ce=(Aloss*dble(dens)**2)/(Eth0*Psc)  ! cgs
              !
              !  apply cooling to primitive and conserved variables
              primit(neqdyn,i,j,k)=primit(neqdyn,i,j,k)*exp(-ce*dt)
              !
              !   u(neqdyn,new)=Ekin0+Eth_new
              u(neqdyn,i,j,k)=u(neqdyn,i,j,k)-Eth0+cv*primit(neqdyn,i,j,k)
              !
           end if
           !
        end do
     end do
  end do
  !
  !--------------------------------------------------------------------
  !
#endif
end subroutine coolingdmc
!======================================================================
