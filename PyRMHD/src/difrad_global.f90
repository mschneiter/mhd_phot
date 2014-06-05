!======================================================================
!   Module to implement diffuse radiative trasnport using ray tracing
!======================================================================
#ifdef RADDIFF_GLOBAL
  module difrad_global
    use parameters
    use globals
    implicit none
    real, parameter   :: a0=6.3e-18
    real, allocatable :: nh0(:,:,:), ph(:,:,:), em(:,:,:), aux(:,:,:)

  contains
    !--------------------------------------------------------------------
    !   initializes random number generation 
    !--------------------------------------------------------------------
    subroutine init_rand
      implicit none
      integer :: rand_size
      integer, allocatable, dimension(:) :: rand_seed
      character (len=10) :: system_time
      real :: rtime
    
      allocate(nh0(1:nxtot,1:nytot,1:nztot) )
      allocate( ph(1:nxtot,1:nytot,1:nztot) )
      allocate( em(1:nxtot,1:nytot,1:nztot) )
      allocate(aux(1:nxtot,1:nytot,1:nztot) )

      call random_seed(size=rand_size)
      allocate(rand_seed(1:rand_size))
      call date_and_time(time=system_time)
      read(system_time,*) rtime
      rand_seed=int(rtime*1000.)
#ifdef MPIP
      rand_seed=rand_seed*rank
#endif
      call random_seed(put=rand_seed)    
      deallocate(rand_seed)

    end subroutine init_rand
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    !   calculates the diffuse fotoionization emissivity
    !--------------------------------------------------------------------
    subroutine emdiff(emax)
      implicit none
      real, intent(out) :: emax
      real :: prim(neq)
      real :: T, emaxp, de, emLym, emHeII
      integer :: i ,j , k, err
      !
      !real :: x,y,z,rad,rstar,vol,emstar
      !
      !   clear  variables
      emaxp=0.
      em (:,:,:)=0.
      aux (:,:,:)=0.

      do i=1,nx
         do j=1,ny
            do k=1,nz       
               call uprim(prim, u(:,i,j,k), T)
               !  electron density
               de=max(u(1,i,j,k)-u(neqdin+1,i,j,k),0.)
               !  Lyman cont. emission
               emLym= de*de*1.57e-13*(1.e4/T)**0.52

               !  He II emission
               emHeII= 2.*0.1*de*de*0.101*8.629e-6*exp(-4.72e5/T)/sqrt(T)

               !  Total emissivity
               aux(coords(0)*nx+i, coords(1)*ny+j, coords(2)*nz+k) = emLym !+emHeII
               !
               emaxp=max(emaxp, aux(coords(0)*nx+i,coords(1)*ny+j,coords(2)*nz+k)  )

            end do
         end do
      end do

#ifdef MPIP
      call mpi_allreduce(emaxp, emax, 1, mpi_real_kind, mpi_max, comm2d,err)
      call mpi_allreduce(aux, em, nxtot*nytot* nztot,mpi_real_kind, mpi_sum, comm2d, err)
#else
      em(:,:,:)=emp(:,:,:)
      emax=aux
#endif

      return
    end subroutine emdiff
    !--------------------------------------------------------------------
    
    !--------------------------------------------------------------------
    !   photon random directions
    !----------------_----------------------------------------------------
    subroutine random_direction(xd,yd,zd)
      implicit none
      real, intent (out) :: xd, yd, zd
      real :: r2, w
      real :: ran(3)

      r2=2.
      do while (r2.gt.1.)
         call random_number(ran)
         xd=2.*(ran(1)-0.5)
         yd=2.*(ran(2)-0.5)
         zd=2.*(ran(3)-0.5)
         r2=xd**2+yd**2+zd**2
      end do
      w=1./sqrt(r2)
      xd=xd*w
      yd=yd*w
      zd=zd*w
      !
    end subroutine random_direction
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    !   place photon packets at a "star" surface
    !--------------------------------------------------------------------
    subroutine starsource(Srad,is,js,ks,i,j,k,xd,yd,zd)
      implicit none
      integer, intent (in)  :: Srad, is, js, ks
      integer, intent (out) :: i,  j,  k
      real, intent (out) :: xd, yd, zd
      real :: xs, ys,zs
      real :: rd, rs, w, prod
      real :: ran(6)

      call random_number(ran)
         
      !   get unit vector from star center
      xs=2.*(ran(1)-0.5)
      ys=2.*(ran(2)-0.5)
      zs=2.*(ran(3)-0.5)
      rs=xs**2+ys**2+zs**2
      w=float(Srad)/sqrt(rs)
      xs=xs*w
      ys=ys*w
      zs=zs*w
      !
      !  get the integer position on the source
      i=is+int( xs )
      j=js+int( ys )
      k=ks+int( zs )
      !
      !  get the direction of the ray
      xd=2.*(ran(4)-0.5)
      yd=2.*(ran(5)-0.5)
      zd=2.*(ran(6)-0.5)
      rd=xd**2+yd**2+zd**2
      w=1./sqrt(rd)
      xd=xd*w
      yd=yd*w
      zd=zd*w
      !
      !  reflect rays that point to the star
      prod=( float(i-is)*xd+ float(j-js)*yd + float(k-ks)*zd )
      if (prod < 0. ) then
         xd=-xd
         yd=-yd
         zd=-zd
      endif
      !
    end subroutine starsource
    !--------------------------------------------------------------------
    

    !--------------------------------------------------------------------
    !   photon trajectories
    !   launches a photon from cell (ic,jc,kc) in the (x1,y1,z1) 
    !   direction, with an ionizing rate f
    !--------------------------------------------------------------------
    subroutine photons(ic,jc,kc,x1,y1,z1,f)
      implicit none
      real, intent(in) :: x1,y1,z1
      real, intent(inout) :: f
      integer, intent(in) :: ic, jc, kc
      real :: xa, xb, ya, yb, za, zb
      real :: dl, dxl, dyl, dzl, xl0, yl0, zl0, xl, yl, zl, al, dtau
      real :: fmin
      integer :: iter, i, j ,k 
      !
      xa=1.
      xb=float(nxtot)
      ya=1.
      yb=float(nytot)
      za=1.
      zb=float(nztot)
      !
      !     photon trajectory
      !
      fmin=1.e-10*f
      dl=0.5
      dxl=dl*x1
      dyl=dl*y1
      dzl=dl*z1
      if(ic.ge.1) then
         xl0=float(ic)
         yl0=float(jc)
         zl0=float(kc)
      else
         if(x1.gt.0.) then
            xl0=1.
            al=(xl0-float(ic))/x1
            yl0=al*y1+float(jc)
            zl0=al*z1+float(kc)
         else
            xl0=0.
            yl0=0.
            zl0=0.
         end if
      end if
      !
      xl=xl0+dxl
      yl=yl0+dyl
      zl=zl0+dzl
      !
      do while( f.ge.fmin                  &
#ifndef PERIODX     
              .and.xl.le.xb.and.xl.ge.xa  &
#endif
#ifndef PERIODY
              .and.yl.le.yb.and.yl.ge.ya  &
#endif
#ifndef PERIODZ
              .and.zl.le.zb.and.zl.ge.za  &
#endif
            )
         !
        i=int(xl+0.5)
        j=int(yl+0.5)
        k=int(zl+0.5)
#ifdef PERIODX
      if (i .lt. 1) then
        i=nxtot
        xl=float(nxtot)+dxl
      else if(i.gt.nxtot) then
        i=1
        xl=1.+dxl
      endif 
#else
       xl=xl+dxl
#endif

#ifdef PERIODY
      if (j .lt. 1) then
        print*,'period',j,yl,dyl
        j=nytot
        yl=float(nytot)+dyl
      else if(j.gt.nytot) then
        print*,'period',j,yl,dyl
        j=1
        yl=1.+dyl
      endif 
#else
       yl=yl+dyl
#endif

#ifdef PERIODZ
      if (k .lt. 1) then
        print*,'period',k,zl,dzl
        k=nztot
        zl=float(nztot)+dzl
      else if(k.gt.nztot) then
        k=1
        zl=1.+dzl
      endif 
#else
       zl=zl+dzl
#endif
         !
         dtau=nh0(i,j,k)*a0*dl*dx*rsc
         aux(i,j,k)=aux(i,j,k)+f*(1.-exp(-dtau) )/nh0(i,j,k)
         f=f*exp(-dtau)
        !
      end do
      !
      !
      return
    end subroutine photons
 !--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !  progress bar (takes a number between 1 and tot)
 !--------------------------------------------------------------------
 subroutine progress(j,tot)
   implicit none
   integer(kind=4)::j,k
   integer(kind=4), intent(in) :: tot
   character(len=57)::bar="???% |                                                  |"
   open (unit=6, carriagecontrol='fortran')
   write(unit=bar(1:3),fmt="(i3)") 100*j/tot
   bar(7:56)="."
   do k=1, 50*j/tot
      bar(6+k:6+k)="="
   enddo
   ! print the progress bar.
   write(unit=6,fmt="(a1,a1,a57)") '+',char(13), bar
   return
 end subroutine progress
 !--------------------------------------------------------------------

 !--------------------------------------------------------------------
 !   Upper level wrapper to compute the diffuse photoionization rate
 !--------------------------------------------------------------------
 subroutine diffuse_rad()
   implicit none
   real :: nmax,emax, f, dirx,diry,dirz, radius
   integer :: i,j,k, nrays, err, niter, ic,jc,kc, Srad

   nmax = nxtot*nytot*nztot/10000/np
   !resets the photoionization rate, nh0 and, aux
   ph (:,:,:)=0.
   aux(:,:,:)=0.
   nh0(:,:,:)=0.
   !    
   !   retrieve the density
   aux( coords(0)*nx+1:(coords(0)*nx+nx), &
        coords(1)*ny+1:(coords(1)*ny+ny), &
        coords(2)*nz+1:(coords(2)*nz+nz) ) = u (neqdin+1,1:nx,1:ny,1:nz)

   call mpi_allreduce (aux, nh0, nxtot*nytot*nztot, &
        mpi_real_kind, mpi_sum, comm2d, err)

   !   computes the emissivity at each cell
   !call emdiff(emax)

   aux (:,:,:)=0.

   !   fire the photon torpedoes!
   
   nrays=100000
   !   posicion de la fuente ionizante
   ic=nxtot/2
   jc=nytot/2
   kc=nztot/2

   !   impose the photoionizing field of a star
   f=1.173e12*a0/float(np)/float(nrays)   !  flujo de fotones ionizantes
!   Srad=int(Rsw/dx)   ! needs to load star module
   srad=7
   do niter=1, nrays 
      ! gets the position on the star sirface and direction of
      ! the photon packet to be traced
      call starsource(Srad,ic,jc,kc,i,j,k,dirx,diry,dirz)
      ! trace the photon
      call photons(i,j,k,dirx,diry,dirz,f)
      !if (rank.eq.0) call progress(niter,nrays)      
   end do

   
!!$   The following is fo the diffuse radiation
!!$   do i=1,nxtot
!!$      if (rank.eq.0) call progress(i,nxtot)
!!$      do j=1,nytot
!!$         do k=1,nztot
!!$            !if (em(i,j,k).gt.(1e-10)) then
!!$            !   nrays=nmax*int(10.*em(i,j,k)/emax)
!!$            !   nrays=max(nrays,5)
!!$            !
!!$            !   do niter=1,nrays
!!$            !     f=em(i,j,k)/float(nrays*np)
!!$            !     call random_direction(dirx,diry,dirz)
!!$            !     call photons(i,j,k,dirx,diry,dirz,f) 
!!$            !   end do
!!$            !
!!$            !end if
!!$            
!!$            !   here i impose the plane parallel ionizing flux
!!$            !if (i.eq.10) then
!!$            !   f=9.E8*a0/float(np)
!!$            !   call random_direction(dirx,diry,dirz)
!!$            !   call photons(i,j,k,dirx,diry,dirz,f) 
!!$            !   call photons(i,j,k,1.,0.,0.,f)
!!$            !end if
!!$
!!$            
!!$            radius=sqrt( float( i-ic)**2 + float(j-jc)**2 + float(k-kc)**2)
!!$            if (radius .le.0.5 ) then
!!$               f=1.173e12*a0/float(np)   !  flujo de fotones ionizantes
!!$               !
!!$               do niter=1, nrays 
!!$                  !   obtiene la direccion al azar
!!$                  call random_direction(dirx,diry,dirz)
!!$                  !  hace el ray tracing de un "foton" (i.e. photon packet) 
!!$                  !   inyectado en i,j,k en la direccion 
!!$                  !  dirx, diry, dirz, con un lfijo de fotones ionizantes f
!!$                  call photons(i,j,k,dirx,diry,dirz,f) 
!!$
!!$               end do
!!$               
!!$            endif
!!$
!!$         end do
!!$      end do
!!$   end do

#ifdef MPIP   
   call mpi_allreduce (aux, ph, nxtot*nytot*nztot, &
        mpi_real_kind, mpi_sum, comm2d, err)
#else
   ph(:,:,:)=aux(:,:,:)
#endif
  !if (rank.eq.0) then
    !print*,ph(1:3,1,1)
    !print*,ph(1:3,1,2)
    !print*,ph(1:3,1,3)
  !endif
   return
 end subroutine diffuse_rad
 !--------------------------------------------------------------------
end module difrad_global
#endif
