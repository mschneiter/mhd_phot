!=========================================================================
!   calculates HLLD fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses primit for the 1st half of timestep (first order)
!          --> 2, uses primit for second order timestep
!=========================================================================
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
    real, dimension(neq)               :: pp
    real, dimension(neq)               ::fL, fR, uL, uR
    real, dimension(neq)               :: uu, uuk, ust
    real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sM
    real :: pTL, pTR, Bx, signBx
    real :: slmsM, srmsM, rhostl, rhostr, sstl, sstr
    real :: pst, el, er, denl, denr, sMmul, sMmur
    real :: vstl, wstl, bystl, bzstl, estl, vdotbl, vstdotbstl
    real :: vstr, wstr, bystr, bzstr, estr, vdotbr, vstdotbstr
    real :: sMmsstl, sMmsstr
    real :: dd, vstst, wstst, bystst, bzstst
    real ::  vststdotbstst, eststl, eststr
    integer :: err
    !
    call cfastX(priml,csl)
    call cfastX(primr,csr)
    !
    sr=max(priml(2)+csl,primr(2)+csr)
    sl=min(priml(2)-csl,primr(2)-csr)
    !
    !=====================================================================
    ! UL region 
    !=====================================================================
    if (sl > 0) then
       call primf(priml,ff)
       return
    endif
    !
    !=====================================================================
    ! UR region 
    !=====================================================================
    if (sr < 0) then
       call primf(primr,ff)
       return
    endif
    ! 
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    Bx= 0.5* (primL(6)+primR(6) )
    signBx= sign(1.,Bx)

    !  Total pressure
    pTL=primL(5) + 0.5*( bx**2+primL(7)**2+primL(8)**2 )
    pTR=primR(5) + 0.5*( bx**2+primR(7)**2+primR(8)**2 )


    slmul=sl-priml(2)
    srmur=sr-primr(2)

    rholul=priml(1)*priml(2)
    rhorur=primr(1)*primr(2)

    sM = (srmur*rhorur-slmul*rholul-pTR+pTL)/( srmur*primr(1)-slmul*priml(1) )
    !
    srmsM=sr-sM
    slmsM=sl-sM
    !
    rhostl=priml(1)*slmul/slmsM      !rhoL*
    rhostr=primr(1)*srmur/srmsM      !rhoR*
    !
    sstl=sM - abs(bx)/sqrt(rhostl)  !SL*
    sstr=sM + abs(bx)/sqrt(rhostr)  !SR*

    pst= (srmur*primr(1)*pTL - slmul*priml(1)*pTR               & 
          + priml(1)*primr(1)*srmur*slmul*(primr(2)-priml(2)) ) &
          /( srmur*primr(1)-slmul*priml(1) )
    !
    !=====================================================================
    ! UL* region 
    !=====================================================================
    !
    if(sstl >= 0) then
      !
      el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
                 +0.5*(bx**2+priml(7)**2+priml(8)**2)
      !
      sMmuL=sM - priml(2)
      denl=priml(1)*slmul*slmsM-bx**2
      !
      if(denl == 0) then
        call primu(primL,uuk)
        stop
      else

        vstl = priml(3) - bx*priml(7)*sMmul/denl                      !vL*
        wstl = priml(4) - bx*priml(8)*sMmul/denl                      !wL*
        bystl= priml(7)*( priml(1)*slmul**2 - bx**2 )/denl            !byL*
        bzstl= priml(8)*( priml(1)*slmul**2 - bx**2 )/denl            !bzL*
        !
        vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
        vstdotbstl= sM*bx + vstl*bystl + wstl*bzstl                     !vL* dot BL*  
        !
        estl= ( slmul*el -pTL*priml(2) +pst*sM +bx*(vdotbl-vstdotbstl) )/slmsM !eL*
        !
        ff(1) = rhostl*sM
        ff(2) = rhostl*SM**2+pst-bx**2
        ff(3) = rhostl*sM*vstl-bx*bystl
        ff(4) = rhostl*sM*wstl-bx*bzstl
        ff(5) = sM*(estl+pst)-bx*(vstdotbstl)
        ff(6) = 0.
        ff(7) = bystl*sM-bx*vstl
        ff(8) = bzstl*sM-bx*wstl
        !
#ifdef PASSIVES
        ff(neqdyn+1:neq)=sM*priml(neqdyn+1:neq)*slmul/slmsM
#endif
      !
      endif
      !
    return
    endif
    !
    !=====================================================================
    ! UR* region 
    !=====================================================================
    !
    if(sstr <= 0) then
      !
      er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
                 +0.5*(bx**2+primr(7)**2+primr(8)**2)
      !
      sMmuR=sM - primr(2)
      denr=primr(1)*srmur*sRmsM-bx**2
      !
      if(denr == 0) then
        call primu(primR,uuk)
        stop
      else

        vstr = primr(3) - bx*primr(7)*sMmuR/denr                      !vR*
        wstr = primr(4) - bx*primr(8)*sMmuR/denr                      !wR*
        bystr= primr(7)*( primr(1)*srmur**2 - bx**2 )/denr            !byR*
        bzstr= primr(8)*( primr(1)*srmur**2 - bx**2 )/denr            !bzR*
        !
        vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
        vstdotbstr= sM*bx + vstr*bystr + wstr*bzstr                     !vR* dot BR*  
        !
        estr= ( srmur*er -pTR*primr(2) +pst*sM +bx*(vdotbr-vstdotbstr) )/srmsM !eR*
        !
        ff(1) = rhostr*sM
        ff(2) = rhostr*SM**2+pst-bx**2
        ff(3) = rhostr*sM*vstr-bx*bystr
        ff(4) = rhostr*sM*wstr-bx*bzstr
        ff(5) = sM*(estr+pst)-bx*(vstdotbstr)
        ff(6) = 0.
        ff(7) = bystr*sM-bx*vstr
        ff(8) = bzstr*sM-bx*wstr
        !
#ifdef PASSIVES
        ff(neqdyn+1:neq)=sM*primr(neqdyn+1:neq)*srmur/srmsM
#endif
    endif
      !
    return
    endif
    !
    !=====================================================================  
    !
    !   All this are needed on both the UL** and UR** regions
    sMmul= sM - priml(2)
    sMmur= sM - primr(2)

    denl=priml(1)*slmul*slmsM-bx**2
    denr=primr(1)*srmur*srmsM-bx**2
    !
    if(denl == 0) then
      vstl =priml(3)
      wstl =priml(4)
      bystl=priml(7)
      bzstl=priml(8)
      stop
    else
      vstl = priml(3) - bx*priml(7)*sMmul/denl                      !vL*
      wstl = priml(4) - bx*priml(8)*sMmul/denl                      !wL*
      bystl= priml(7)*( priml(1)*slmul**2 - bx**2 )/denl            !byL*
      bzstl= priml(8)*( priml(1)*slmul**2 - bx**2 )/denl            !bzL*
    endif

    if(denr == 0) then
      vstr =primr(3)
      wstr =primr(4)
      bystr=primr(7)
      bzstr=primr(8)
      stop
      else
      vstr = primr(3) - bx*primr(7)*sMmuR/denr                      !vR*
      wstr = primr(4) - bx*primr(8)*sMmuR/denr                      !wR*
      bystr= primr(7)*( primr(1)*srmur**2 - bx**2 )/denr            !byR*
      bzstr= primr(8)*( primr(1)*srmur**2 - bx**2 )/denr            !bzR*
    endif

    dd=sqrt(rhostl)+sqrt(rhostr)

    vstst =(sqrt(rhostl)*vstl + sqrt(rhostr)*vstr + (bystr-bystl)*signBx )/dd  !v** 
    wstst =(sqrt(rhostl)*wstl + sqrt(rhostr)*wstr + (bzstr-bzstl)*signBx )/dd  !w**

    bystst=(sqrt(rhostl)*bystr + sqrt(rhostr)*bystl +  &      !by**
            sqrt(rhostl*rhostr)*(vstr-vstl)*signBx )/dd                            
   
    bzstst=(sqrt(rhostl)*bzstr + sqrt(rhostr)*bzstl +  &      !bz**
            sqrt(rhostl*rhostr)*(wstr-wstl)*signBx )/dd

    vststdotbstst= sM*bx + vstst*bystst + wstst*bzstst        !v** dot B**
    !
    !=====================================================================
    ! UL** region 
    !=====================================================================
    !
    if(sM >= 0 ) then
      !
      el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
                  +0.5*(bx**2+priml(7)**2+priml(8)**2)
      !
      vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
      vstdotbstl= sM*bx+ vstl*bystl + wstl*bzstl                      !vL* dot BL*  
      !
      estl= ( slmul*el -pTL*priml(2) +pst*sM +bx*(vdotbl-vstdotbstl) )/slmsM !eL*
      !
      eststl= estl - sqrt(rhostl)*(vstdotbstl-vststdotbstst)*signBx !eL**
      ! 
      ff(1) = rhostl*sM
      ff(2) = rhostl*SM**2+pst-bx**2
      ff(3) = rhostl*sM*vstst-bx*bystst
      ff(4) = rhostl*sM*wstst-bx*bzstst
      ff(5) = sm*(eststl+pst)-bx*(vststdotbstst)
      ff(6) = 0.
      ff(7) = bystst*sM-bx*vstst
      ff(8) = bzstst*sM-bx*wstst
      !
#ifdef PASSIVES
      ff(neqdyn+1:neq) = sM*priml(neqdyn+1:neq)*slmul/slmsM
#endif
    !
    return
    endif
    !
    !=====================================================================
    ! UR** region 
    !=====================================================================
    !
    if(sM <= 0 ) then
      !
      er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
                  +0.5*(bx**2+primr(7)**2+primr(8)**2)
      !
      vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
      vstdotbstr= sM*bx+ vstr*bystr + wstr*bzstr                      !vR* dot BR*  
      !
      estr= ( srmur*er -pTR*primr(2) +pst*sM +bx*(vdotbr-vstdotbstr) )/srmsM !eR*
      !
      eststr= estr + sqrt(rhostr)*(vstdotbstr-vststdotbstst)*signBx !eR**
      !
      ff(1) = rhostr*sM
      ff(2) = rhostr*SM**2+pst-bx**2
      ff(3) = rhostr*sM*vstst-bx*bystst
      ff(4) = rhostr*sM*wstst-bx*bzstst
      ff(5) = sm*(eststr+pst)-bx*(vststdotbstst)
      ff(6) = 0.
      ff(7) = bystst*sM-bx*vstst
      ff(8) = bzstst*sM-bx*wstst
      !
#ifdef PASSIVES
      ff(neqdyn+1:neq) = sM*primr(neqdyn+1:neq)*srmur/srmsM)
#endif
      !
    return
    endif
    !   
    !=====================================================================
    !
    print'(a,5es12.3)', 'Error in HLLD routine', sM,sl,sr, csl, csr
    stop
    !
    end subroutine primfhlld
#endif
end subroutine hlldfluxes
  !=======================================================================
