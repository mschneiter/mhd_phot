!======================================================================
!   Cooling routine with tabulated values from CHIANTI
!   The cooling is turned of bellow 10^4 K
!   And the primitives are updated in this routine as well
!======================================================================
subroutine coolingchi(dt)
#ifdef COOLINGCHI
  use parameters
  use globals
  use chi_module
  implicit none
  real,    intent(in)  :: dt
  real                 :: T ,Eth0, dens
  real, parameter      :: Tmin=10000.
  real (kind=8)        :: ALOSS, Ce
  integer              :: i, j, k
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
              Aloss=coolchi(T)
              dens=primit(1,i,j,k)
              Ce=(Aloss*dble(dens)**2)/(Eth0*Psc)  ! cgs
              !
              !  apply cooling to primitive and conserved variables
              primit(neqdyn,i,j,k)=primit(neqdyn,i,j,k)*exp(-ce*dt)
              !
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
end subroutine coolingchi
!======================================================================
