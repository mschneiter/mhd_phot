!=======================================================================
!   calculates HLLD fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses primit for the 1st half of timestep (first order)
!          --> 2, uses primit for second order timestep
!=======================================================================
subroutine hlldfluxes(choice)
#ifdef HLLD
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
              call primfhlld(priml,primr,ff)
              f(:,i,j,k)=ff(:)
              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, jp,k )
              call swapy(priml,neq)          !swaps primL for L state
              call swapy(primr,neq)          !swaps primR for R state 
              !
              call primfhlld(priml,primr,ff)  !gets fluxes (swapped)
              call swapy(ff,neq)             !swaps back the fluxes
              g(:,i,j,k)=ff(:)
              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, j, kp)
              call swapz(priml,neq)
              call swapz(primr,neq)
              !
              call primfhlld(priml,primr,ff)
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
              call primfhlld(priml,primr,ff)
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
              call primfhlld(priml,primr,ff)
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
              call primfhlld(priml,primr,ff)
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
  !   calculates the F HLLD fluxes from the primitive variables
  !=======================================================================
  subroutine primfhlld(priml,primr,ff)
    use parameters
    use sound 
    implicit none
    real, dimension(neq),intent(in   ) :: priml, primr   
    real, dimension(neq),intent(inout) :: ff
    real, dimension(neq)               :: fL, fR, uL, uR
    real, dimension(neq)               :: uu, uuk, ust
    real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sst
    real :: slmsst, srmsst, rhostl, rhostr, sstl, sstr
    real :: pst, el, er, denl, denr, sstmul, sstmur
    real :: vstl, wstl, bystl, bzstl, estl, vdotbl, vstdotbstl
    real :: vstr, wstr, bystr, bzstr, estr, vdotbr, vstdotbstr
    real :: sstmsstl, sstmsstr
    real :: dd, vstst, wstst, bystst, bzstst
    real ::  vststdotbstst, eststl, eststr
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
    slmul=sl-priml(2) 
    srmur=sr-primr(2)
    rholul=priml(1)*priml(2)
    rhorur=primr(1)*primr(2)
    !
    sst = (srmur*rhorur-slmul*rholul-primr(5)+priml(5) )        &                           !SM
         / (srmur*primr(1)-slmul*priml(1) )
    !
    !
    sstl=sst-abs(priml(6))/sqrt(rhostl)                                                     !SL*
    sstr=sst+abs(primr(6))/sqrt(rhostr)                                                     !SR*
    !
    pst= (srmur*primr(1)*priml(5) - slmul*priml(1)*primr(5)     &                           !p* 
          + priml(1)*primr(1)*srmur*slmul*(primr(2)-priml(2)) ) &
          /(srmur*primr(1)-slmul*priml(1) )
    !
    !--------------------------------------------------------------------------
    !Variables for the left star region   
    !--------------------------------------------------------------------------
    !
    slmsst=sl-sst 
    rhostl=priml(1)*slmul/slmsst                                                            !rhoL*
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &                     !eL
                 +0.5*(priml(6)**2+priml(7)**2+priml(8)**2)
    !
    sstmul=sst - priml(2)
    denl=priml(1)*slmul*slmsst-priml(6)**2
    !
    vstl = priml(3) - priml(6)*priml(7)*sstmul/denl                                         !vL*
    wstl = priml(4) - priml(6)*priml(8)*sstmul/denl                                         !wL*
    bystl= (priml(7)*priml(1)*slmul**2-priml(6)**2)/denl                                    !byL*
    bzstl= (priml(8)*priml(1)*slmul**2-priml(6)**2)/denl                                    !bzL*
    !
    vdotbl=priml(2)*priml(6)+priml(3)*priml(7)+priml(4)*priml(8)                            !vL dot BL
    vstdotbstl=priml(2)*priml(6)+ vstl*bystl+wstl*bzstl                                     !vL* dot BL*  
    !
    estl= slmul*el-priml(5)*priml(2)+pst*sst+priml(6)*(vdotbl-vstdotbstl)/slmsst            !eL*
    !
    !-------------------------------------------------------------------------------
    !Variables for the right star region
    !------------------------------------------------------------------------------
    !
    srmsst=sr-sst 
    rhostr=primr(1)*srmur/srmsst                                                            !rhoR*
    !
    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &                     !eR
                 +0.5*(primr(6)**2+primr(7)**2+primr(8)**2)
    !
    sstmur=sst - primr(2)
    denr=primr(1)*srmur*srmsst-primr(6)**2
    !
    vstr = primr(3) - primr(6)*primr(7)*sstmur/denr                                         !vR*
    wstr = primr(4) - primr(6)*primr(8)*sstmur/denr                                         !wR*
    bystr= (primr(7)*primr(1)*srmur**2-primr(6)**2)/denr                                    !byR*
    bzstr= (primr(8)*primr(1)*srmur**2-primr(6)**2)/denr                                    !bzR*
    !
    vdotbr=primr(2)*primr(6)+primr(3)*primr(7)+primr(4)*primr(8)                            !vR dot BR
    vstdotbstr=primr(2)*primr(6)+ vstr*bystr+wstr*bzstr                                     !vR* dot BR*
    !
    estr= srmur*er-primr(5)*primr(2)+pst*sst+primr(6)*(vdotbr-vstdotbstr)/srmsst            !eR*
    !
    !
    sstmsstl=sst-sstl                                                                       !s*-sL*
    sstmsstr=sst-sstr                                                                       !s*-sR*
    !-------------------------------------------------------------------------------
    !Variables for the star star region
    !------------------------------------------------------------------------------
    !
    dd=sqrt(rhostl)+sqrt(rhostr)
    vstst =(sqrt(rhostl)*vstl + sqrt(rhostr)*vstr + (bystr-bystl)*sign(1.0,priml(6)))/dd    !v** 
    wstst =(sqrt(rhostl)*wstl + sqrt(rhostr)*wstr + (bzstr-bzstl)*sign(1.0,priml(6)))/dd    !w**
    bystst=(sqrt(rhostl)*bystr + sqrt(rhostr)*bystl +  &                                    !by**
            sqrt(rhostl*rhostr)*(vstr-vstl)*sign(1.0,priml(6)))/dd                            
    bzstst=(sqrt(rhostl)*bzstr + sqrt(rhostr)*bzstl +  &                                    !bz**
            sqrt(rhostl*rhostr)*(wstr-wstl)*sign(1.0,primr(6)))/dd
    ! 
    !---------------------------------------------------------------------------
    !============================================
    ! left star star region 
    !============================================
    if(sstmsstl .ge. 0 ) then
    ! 
       vststdotbstst= sst*priml(6)+vstst*bystst + wstst*bzstst                              !v** dot B**
       eststl= estl - sqrt(rhostl)*(vstdotbstl-vststdotbstst)*sign(1.0,priml(6))            !eL**
    !  
       uuk(1)=rhostl
       uuk(2)=rhostl*sstl
       uuk(3)=rhostl*vstst
       uuk(4)=rhostl*wstst
       uuk(5)=eststl
       uuk(6)=priml(6)
       uuk(7)=bystst 
       uuk(8)=bzstst   
    !      
       ust(1)=rhostl
       ust(2)=rhostl*sstl
       ust(3)=rhostl*vstl
       ust(4)=rhostl*wstl
       ust(5)=estl
       ust(6)=priml(6) 
       ust(7)=bystl     
       ust(8)=bzstl   
    !
#ifdef PASSIVES
       uuk(neqdyn+1:neq)=priml(neqdyn+1:neq)*slmul/slmsst
       ust(neqdyn+1:neq)=priml(neqdyn+1:neq)*slmul/slmsst
#endif
    !
       call primf(priml,fL)
       call primu(priml,uL)
    !
       ff(:)=fL(:)+sstl*uuk(:)-(sstl-sl)*ust(:)-sl*uL(:)
    !    
    return
    endif
    !-----------------------------------------------------------------------
    !============================================
    ! right star star region 
    !============================================
    if(sstmsstr .le. 0 ) then
    !
       vststdotbstst= sst*primr(6)+vstst*bystst + wstst*bzstst                              !v** dot B**
       eststr= estr + sqrt(rhostr)*(vstdotbstr-vststdotbstst)*sign(1.0,primr(6))            !eL**
    !
       uuk(1)=rhostr
       uuk(2)=rhostr*sst
       uuk(3)=rhostr*vstst
       uuk(4)=rhostr*wstst
       uuk(5)=eststr
       uuk(6)=primr(6)
       uuk(7)=bystst 
       uuk(8)=bzstst
    !
       ust(1)=rhostr
       ust(2)=rhostr*sst
       ust(3)=rhostr*vstr
       ust(4)=rhostr*wstr
       ust(5)=estr
       ust(6)=primr(6)    
       ust(7)=bystr    
       ust(8)=bzstr   
    !
#ifdef PASSIVES
       uuk(neqdyn+1:neq)=primr(neqdyn+1:neq)*srmur/srmsst
       ust(neqdyn+1:neq)=primr(neqdyn+1:neq)*srmur/srmsst
#endif
     
       call primf(primr,fR)
       call primu(primr,uR)
    !
       ff(:)=fR(:)+sstr*uuk(:)-(sstr-sr)*ust(:)-sr*uR(:)
    !
        return
    endif
    !-----------------------------------------------------------------------

    !============================================
    ! left star region 
    !============================================
    if(sstl .ge. 0) then
    !
       uuk(1)=rhostl
       uuk(2)=rhostl*sst
       uuk(3)=rhostl*vstl
       uuk(4)=rhostl*wstl
       uuk(5)=estl
       uuk(6)=priml(6)    
       uuk(7)=bystl    
       uuk(8)=bzstl   
    !
#ifdef PASSIVES
       uuk(neqdyn+1:neq)=priml(neqdyn+1:neq)*slmul/slmsst
#endif
 
       call primf(priml,fL)
       call primu(priml,uL)
    !
       ff(:)=fL(:)+sl*uuk(:)-sl*uL(:)
    !
    return
    endif
    !--------------------------------

    !============================================
    ! right star region 
    !============================================
    if(sstr .ge. 0) then
    !
       uuk(1)=rhostr
       uuk(2)=rhostr*sst
       uuk(3)=rhostr*vstr
       uuk(4)=rhostr*wstr
       uuk(5)=estr
       uuk(6)=primr(6)    
       uuk(7)=bystr    
       uuk(8)=bzstr   
    !
#ifdef PASSIVES
       uuk(neqdyn+1:neq)=primr(neqdyn+1:neq)*srmur/srmsst
#endif
    !
       call primf(primr,fR)
       call primu(primr,uR)
    !
       ff(:)=fR(:)+sr*uuk(:)-sr*uR(:)
    !
    return
    endif
    !-----------------------------------------------------------------------
    end subroutine primfhlld
#endif
end subroutine hlldfluxes
  !=======================================================================
