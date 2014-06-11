!=======================================================================
!   calculates HLLE fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses primit for the 1st half of timestep (first order)
!          --> 2, uses primit for second order timestep
!=======================================================================
subroutine hllEfluxes(choice)
#ifdef HLLE
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
              call primfhlle(priml,primr,ff)
              f(:,i,j,k)=ff(:)
              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, jp,k )
              call swapy(priml,neq)          !swaps primL for L state
              call swapy(primr,neq)          !swaps primR for R state 
              !
              call primfhlle(priml,primr,ff)  !gets fluxes (swapped)
              call swapy(ff,neq)             !swaps back the fluxes
              g(:,i,j,k)=ff(:)
              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, j, kp)
              call swapz(priml,neq)
              call swapz(primr,neq)
              !
              call primfhlle(priml,primr,ff)
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
              call primfhlle(priml,primr,ff)
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
              call primfhlle(priml,primr,ff)
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
              call primfhlle(priml,primr,ff)
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
  !   calculates the Fast magnetosonic wave in the X-direcion
  !=======================================================================
 subroutine cfastX(prim,cfX)
  implicit none
  real, intent(in) :: prim(neq)
  real, intent(out) ::cfX
  real :: b2, cs2va2
  !
  b2=prim(6)**2+prim(7)**2+prim(8)**2
  cs2va2 = (gamma*prim(5)+b2)/prim(1)   ! cs^2 + ca^2

  cfx=sqrt(0.5*(cs2va2+sqrt(cs2va2**2-4.*gamma*prim(5)*prim(6)**2/prim(1)/prim(1) ) ) )

  !  return
end subroutine cfastX
  !
  !=======================================================================
  !   calculates the F HLLE fluxes from the primitive variables
  !=======================================================================
  subroutine primfhlle(priml,primr,ff)
    use parameters
    use sound 
    implicit none
    real, dimension(neq),intent(in   ) :: priml, primr
    real, dimension(neq),intent(inout) :: ff
    real, dimension(neq)               :: uR, uL, fL, fR
    real :: csl, csr, sl, sr
    !
    call cfastX(priml,csl)
    call cfastX(primr,csr)
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
    end subroutine primfhlle
#endif
end subroutine hllefluxes
  !=======================================================================
