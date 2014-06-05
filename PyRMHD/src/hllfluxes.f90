!=======================================================================
!   calculates HLL fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses primit for the 1st half of timestep (first order)
!          --> 2, uses primit for second order timestep
!=======================================================================
subroutine hllfluxes(choice)
#ifdef HLL
  use parameters
  use globals
  implicit none
  integer, intent(in) :: choice
  integer :: i, j, k, ip, jp, kp, im, jm, km, ip2, jp2, kp2
  real, dimension(neq) :: priml, primr, primll, primrr, ff, uu
  !
  select case(choice)
     !------------------------------------------------------------------
  case(1)        ! 1st half timestep
     !
     do i=0,nx
        do j=0,ny
           do k=0,nz
              !
              ip=i+1
              jp=j+1
              kp=k+1
              !
              !------- x direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,ip,j ,k )
              !
              call primfhll(priml,primr,ff)  !!!!!!!REVISAR!!!!!!!!!!!!!
              f(:,i,j,k)=ff(:)
              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, jp,k )
              call swapy(priml,neq)          !SWAPEA PRIML PARA ESTADO L   
              call swapy(primr,neq)          !SWAPEA PRIMR PARA ESTADO R 
              !
              call primfhll(priml,primr,ff) !CALCULA FLUJOS SWAPEADOS
              call swapy(ff,neq)             !ACOMODA FLUJOS 
              g(:,i,j,k)=ff(:)
              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, j, kp)
              call swapz(priml,neq)
              call swapz(primr,neq)
              !
              call primfhll(priml,primr,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)
              !
           end do
        end do
     end do
     !------------------------------------------------------------------
  case (2)   !  2nd half timestep
     !------------------------------------------------------------------
     !
     do i=0,nx
        do j=0,ny
           do k=0,nz
              !
              ip=i+1
              ip2=i+2
              im=i-1
              jp=j+1
              jp2=j+2
              jm=j-1
              kp=k+1
              kp2=k+2
              km=k-1
              !
              !------- x direction ------------------------------------
              priml (:)=primit(:,i,  j,k )
              primr (:)=primit(:,ip, j,k )
              primll(:)=primit(:,im, j,k )
              primrr(:)=primit(:,ip2,j,k )
              call limiter(primll,priml,primr,primrr,neq)
              !
              call primfhll(priml,primr,ff)
              f(:,i,j,k)=ff(:)
              !------- y direction ------------------------------------
              priml (:)=primit(:,i,j  ,k )
              primr (:)=primit(:,i,jp ,k )
              primll(:)=primit(:,i,jm ,k )
              primrr(:)=primit(:,i,jp2,k )
              call swapy(priml,neq)
              call swapy(primr,neq)
              call swapy(primll,neq)
              call swapy(primrr,neq)
              call limiter(primll,priml,primr,primrr,neq)
              !
              call primfhll(priml,primr,ff)
              call swapy(ff,neq)
              g(:,i,j,k)=ff(:)
              !------- z direction ------------------------------------
              priml (:)=primit(:,i,j,k  )
              primr (:)=primit(:,i,j,kp )
              primll(:)=primit(:,i,j,km )
              primrr(:)=primit(:,i,j,kp2)
              call swapz(priml,neq)
              call swapz(primr,neq)
              call swapz(primll,neq)
              call swapz(primrr,neq)
              call limiter(primll,priml,primr,primrr,neq)
              !
              call primfhll(priml,primr,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)
              !
           end do
        end do
     end do
     !----------------------------------------------------------------
  end select
  !
contains
  !
  !=======================================================================
  !   calculates the F HLL fluxes from the primitive variables
  !=======================================================================
  subroutine primfhll(priml,primr,ff)
    use parameters
    use sound 
    implicit none
    real, dimension(neq),intent(in   ) :: priml, primr
    real, dimension(neq),intent(inout) :: ff
    real, dimension(neq)               :: uR, uL, fL, fR
    real :: csl, csr, sl, sr
    !
    call csound(priml(5),priml(1),csl)
    call csound(primr(5),primr(1),csr)
    !
    sr=max(priml(2)+csl,primr(2)+csr)
    sl=min(priml(2)-csl,primr(2)-csr)
    !
    if (sl > 0) then
       call primf(priml,ff)
       return
    endif
    !--------------------------------
    if (sr < 0) then
       call primf(primr,ff)
       return
    endif
    !--------------------------------
    call primf(priml,fL)
    call primf(primr,fR)
    call primu(priml,uL)
    call primu(primr,uR)
    !
    ff(:)=(sr*fL(:)-sl*fR(:)+sl*sr*(uR(:)-uL(:)))/(sr-sl)
    return
    !--------------------------------
    end subroutine primfhll
#endif
end subroutine hllfluxes
  !=======================================================================
