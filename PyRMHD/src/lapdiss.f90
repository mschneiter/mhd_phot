!======================================================================
!   Lapidus viscosity (Colella Woodward 1984)
!   choice=1 applied with u
!======================================================================
subroutine lapdiss(choice)
  use parameters
  real :: vx, vy, vxm, vxp, vym, vyp, vzp, vzm
  integer :: i, j, ip, im, jp, jm, kp, km, choice
  !
  select case(choice)
     !------------------------------------------------------------------
  case(1)        ! uses u 1st half timestep
   do i=0,nx
    do j=0,ny
     do k=0,nz
        !
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        kp=k+1
        km=k-1
        !
        vxp =   u(2,ip,j ,k )/u(1,ip,j ,k )
        vxm =   u(2,i ,j ,k )/u(1,i, j ,k )
        vyp = ( u(3,i ,jp,k )/u(1,i ,jp,k )+u(3,ip,jp,k )/u(1,ip,jp,k ) )/2.
        vym = ( u(3,i ,jm,k )/u(1,i ,jm,k )+u(3,im,jm,k )/u(1,im,jm,k ) )/2.
        vzp = ( u(4,i ,j ,kp)/u(1,i ,j ,kp)+u(4,ip,j ,kp)/u(1,ip,j ,kp) )/2.
        vzm = ( u(4,i ,j ,km)/u(1,i ,j ,km)+u(4,im,j ,km)/u(1,im,j ,km) )/2.
        !
        vx=eta*max(0.,vxm-vxp+(vym-vyp)/2.+(vzm-vzp)/2.)
        !
        vxp = ( u(2,ip,j ,k )/u(1,ip,j ,k )+u(2,ip,jp,k )/u(1,ip,jp,k ) )/2.
        vxm = ( u(2,im,j ,k )/u(1,im,j ,k )+u(2,im,jm,k )/u(1,im,jm,k ) )/2.
        vyp =   u(3,i ,jp,k )/u(1,i ,jp,k )
        vym =   u(3,i ,j ,k )/u(1,i ,j ,k )
        vzp = ( u(4,i ,j ,kp)/u(1,i ,j ,kp)+u(4,i ,jp,kp)/u(1,i ,jp,kp) )/2.
        vzm = ( u(4,i ,j ,km)/u(1,i ,j ,km)+u(4,i ,jm,km)/u(1,i ,jm,km) )/2.
        !
        vy=eta*max(0.,(vxm-vxp)/2.+vym-vyp+(vzm-vzp)/2.)
        !
        vxp = ( u(2,ip,j ,k )/u(1,ip,j ,k )+u(2,ip,j ,kp)/u(1,ip,j ,kp) )/2.
        vxm = ( u(2,im,j ,k )/u(1,im,j ,k )+u(2,im,j ,km)/u(1,im,j ,km) )/2.
        vyp = ( u(3,i ,jp,k )/u(1,i ,jp,k )+u(3,i ,jp,kp)/u(1,i ,jp,kp) )/2.
        vym = ( u(3,i ,jm,k )/u(1,i ,jm,k )+u(3,i ,jm,km)/u(1,i ,jm,km) )/2.
        vyp =   u(4,i ,j ,kp)/u(1,i ,j ,kp)
        vym =   u(4,i ,j ,k )/u(1,i ,j ,k )
        !
        vz=eta*max(0.,(vxm-vxp)/2.+(vym-vyp)/2.+vzm-vzp)
        !
        f(:,i,j,k)=f(:,i,j,k)+vx*(u(:,i,j,k)-u(:,ip,j,k))
        g(:,i,j,k)=g(:,i,j,k)+vy*(u(:,i,j,k)-u(:,i,jp,k))
        h(:,i,j,k)=h(:,i,j,k)+vz*(u(:,i,j,k)-u(:,i,j,kp))
        !
     end do
    end do
   end do
  !-----------------------------------------------------------------
  case(2)        ! uses up 2nd half timestep
    !
    do i=0,nx
     do j=0,ny
      do k=0,nz
         !
         ip=i+1
         im=i-1
         jp=j+1
         jm=j-1
         kp=k+1
         km=k-1
         !
         vxp =   up(2,ip,j ,k )/up(1,ip,j ,k )
         vxm =   up(2,i ,j ,k )/up(1,i, j ,k )
         vyp = ( up(3,i ,jp,k )/up(1,i ,jp,k )+up(3,ip,jp,k )/up(1,ip,jp,k ) )/2.
         vym = ( up(3,i ,jm,k )/up(1,i ,jm,k )+up(3,im,jm,k )/up(1,im,jm,k ) )/2.
         vzp = ( up(4,i ,j ,kp)/up(1,i ,j ,kp)+up(4,ip,j ,kp)/up(1,ip,j ,kp) )/2.
         vzm = ( up(4,i ,j ,km)/up(1,i ,j ,km)+up(4,im,j ,km)/up(1,im,j ,km) )/2.
         !
         vx=eta*max(0.,vxm-vxp+(vym-vyp)/2.+(vzm-vzp)/2.)
         !
         vxp = ( up(2,ip,j ,k )/up(1,ip,j ,k )+up(2,ip,jp,k )/up(1,ip,jp,k ) )/2.
         vxm = ( up(2,im,j ,k )/up(1,im,j ,k )+up(2,im,jm,k )/up(1,im,jm,k ) )/2.
         vyp =   up(3,i ,jp,k )/up(1,i ,jp,k )
         vym =   up(3,i ,j ,k )/up(1,i ,j ,k )
         vzp = ( up(4,i ,j ,kp)/up(1,i ,j ,kp)+up(4,i ,jp,kp)/up(1,i ,jp,kp) )/2.
         vzm = ( up(4,i ,j ,km)/up(1,i ,j ,km)+up(4,i ,jm,km)/up(1,i ,jm,km) )/2.
         !
         vy=eta*max(0.,(vxm-vxp)/2.+vym-vyp+(vzm-vzp)/2.)
         !
         vxp = ( up(2,ip,j ,k )/up(1,ip,j ,k )+up(2,ip,j ,kp)/up(1,ip,j ,kp) )/2.
         vxm = ( up(2,im,j ,k )/up(1,im,j ,k )+up(2,im,j ,km)/up(1,im,j ,km) )/2.
         vyp = ( up(3,i ,jp,k )/up(1,i ,jp,k )+up(3,i ,jp,kp)/up(1,i ,jp,kp) )/2.
         vym = ( up(3,i ,jm,k )/up(1,i ,jm,k )+up(3,i ,jm,km)/up(1,i ,jm,km) )/2.
         vyp =   up(4,i ,j ,kp)/up(1,i ,j ,kp)
         vym =   up(4,i ,j ,k )/up(1,i ,j ,k )
         !
         vz=eta*max(0.,(vxm-vxp)/2.+(vym-vyp)/2.+vzm-vzp)
         !
         f(:,i,j,k)=f(:,i,j,k)+vx*(up(:,i,j,k)-up(:,ip,j,k))
         g(:,i,j,k)=g(:,i,j,k)+vy*(up(:,i,j,k)-up(:,i,jp,k))
         h(:,i,j,k)=h(:,i,j,k)+vz*(up(:,i,j,k)-up(:,i,j,kp))
         !
      end do
     end do
    end do
    !------------------------------------------------------------------
  end select
  return
end subroutine lapdiss
!=======================================================================
