!=======================================================================
!   calculates HLLC fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses primit for the 1st half of timestep (first order)
!          --> 2, uses primit for second order timestep
!=======================================================================
subroutine hllcfluxes(choice)
#ifdef HLLC
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
              call primfhllc(priml,primr,ff)  !!!!!!!REVISAR!!!!!!!!!!!!!
              f(:,i,j,k)=ff(:)
              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, jp,k )
              call swapy(priml,neq)          !SWAPEA PRIML PARA ESTADO L   
              call swapy(primr,neq)          !SWAPEA PRIMR PARA ESTADO R 
              !
              call primfhllc(priml,primr,ff) !CALCULA FLUJOS SWAPEADOS
              call swapy(ff,neq)             !ACOMODA FLUJOS 
              g(:,i,j,k)=ff(:)
              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, j, kp)
              call swapz(priml,neq)
              call swapz(primr,neq)
              !
              call primfhllc(priml,primr,ff)
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
              call primfhllc(priml,primr,ff)
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
              call primfhllc(priml,primr,ff)
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
              call primfhllc(priml,primr,ff)
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
  !   calculates the F HLLC fluxes from the primitive variables
  !=======================================================================
  subroutine primfhllc(priml,primr,ff)
    use parameters
    use sound
    implicit none
    real, dimension(neq),intent(in   ) :: priml, primr
    real, dimension(neq),intent(inout) :: ff
    real, dimension(neq)               :: uu, uuk
    real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sst
    real :: rhost,ek
    !
    call csound(priml(5),priml(1),csl)
    call csound(primr(5),primr(1),csr)
    !
    sr=max(priml(2)+csl,primr(2)+csr)
    sl=min(priml(2)-csl,primr(2)-csr)
    !
    if (sl.gt.0) then
       call primf(priml,ff)
       return
    endif
    !--------------------------------
    if (sr.lt.0) then
       call primf(primr,ff)
       return
    endif
    !--------------------------------
    slmul=sl-priml(2)
    srmur=sr-primr(2)
    rholul=priml(1)*priml(2)
    rhorur=primr(1)*primr(2)
    !
    sst = (srmur*rhorur-slmul*rholul-primr(5)+priml(5) )        &  
         / (srmur*primr(1)-slmul*priml(1) )
    !--------------------------------
    if (sst.ge.0.) then
       rhost=priml(1)*(slmul)/(sl-sst)
       ek= 0.5*priml(1)*(priml(2)**2.+priml(3)**2.+priml(4)**2.)+cv*priml(5)
       !
       uuk(1)=rhost
       uuk(2)=rhost*sst
       uuk(3)=rhost*priml(3)
       uuk(4)=rhost*priml(4)
       uuk(5)=rhost*( ek/priml(1)+(sst-priml(2))*(sst+priml(5)/(priml(1)*slmul)) )
       !
#ifdef PMHD
      !uuk(5)= 0.
      uuk(6:8)=rhost*priml(6:8)/priml(1)
#endif
#ifdef PASSIVES
      uuk(neqdyn+1:neq)=rhost*priml(neqdyn+1:neq)/priml(1)
#endif
       !
       call primf(priml,ff)
       call primu(priml,uu)
       ff(:)=ff(:) + sl*( uuk(:)-uu(:) )
       return
    endif
    !--------------------------------
    if (sst.le.0.) then
      rhost=primr(1)*(srmur)/(sr-sst)
      ek= 0.5*primr(1)*(primr(2)**2.+primr(3)**2.+primr(4)**2.)+cv*primr(5)
      !
      uuk(1)=rhost
      uuk(2)=rhost*sst
      uuk(3)=rhost*primr(3)
      uuk(4)=rhost*primr(4)
      uuk(5)=rhost*( ek/primr(1)+(sst-primr(2))*(sst+primr(5)/(primr(1)*srmur)) )
      !
#ifdef PMHD
      !uuk(5)= 0.
      uuk(6:8)=rhost*primr(6:8)/primr(1)
#endif
#ifdef PASSIVES
      uuk(neqdyn+1:neq)=rhost*primr(neqdyn+1:neq)/primr(1)
#endif
       !
       call primf(primr,ff)
       call primu(primr,uu)
       ff(:)=ff(:) + sr*( uuk(:)-uu(:) )
       return
    endif
    !--------------------------------
    !
    print*, 'Error in hllc' 
    stop
    !
  end subroutine primfhllc
#endif
end subroutine hllcfluxes
  !=======================================================================
