!======================================================================
!   Module to implement diffuse radiative trasnport using ray tracing
!======================================================================
#ifdef RADDIFF
module difrad
  use parameters
  use globals
  implicit none
  real, parameter :: a0=6.3e-18
  integer, parameter :: nrays=1000000
  real, allocatable :: ph(:,:,:), em(:,:,:)
  real, allocatable :: photL(:,:,:),photR(:,:,:)
  real, allocatable :: photB(:,:,:),photT(:,:,:)
  real, allocatable :: photO(:,:,:),photI(:,:,:)
  integer :: buffersize(6)
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

    allocate( ph(nx,ny,nz) )
    allocate( em(nx,ny,nz) )

    !allocate buffers for transmission of photons
    !  1st index 1:send 2:recv
    !  2nd index 3xposition + 3xdirection + 1for the F
    !  3rd index a large number can be changed if the code exits
    !  on error
    !  ray positions & directions
    allocate (photL(2,7,nrays))
    allocate (photR(2,7,nrays))
    allocate (photB(2,7,nrays))
    allocate (photT(2,7,nrays))
    allocate (photO(2,7,nrays))
    allocate (photI(2,7,nrays))


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

    do i=1,nx
       do j=1,ny
          do k=1,nz       
             call uprim(prim, u(:,i,j,k), T)
             !  electron density
             de=max(u(1,i,j,k)-u(neqdyn+1,i,j,k),0.)
             !  Lyman cont. emission
             emLym= de*de*1.57e-13*(1.e4/T)**0.52

             !  He II emission (turned off)
             emHeII= 0. ! 2.*0.1*de*de*0.101*8.629e-6*exp(-4.72e5/T)/sqrt(T)

             !  Total emissivity
             em(i,j,k) = emLym +emHeII
             !
             emaxp=max(emaxp, emLym+emHeII)

          end do
       end do
    end do

#ifdef MPIP
    call mpi_allreduce(emaxp, emax, 1, mpi_real_kind, mpi_max, comm2d,err)
#else
    emax=emaxp
#endif

    return
  end subroutine emdiff
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  random vector (with unitary magnutude)
  !----------------_----------------------------------------------------
  subroutine random_versor(xd,yd,zd)
    implicit none
    real, intent (out) :: xd, yd, zd
    real :: r2, w
    real :: ran(3)

    !r2=2.
    !do while (r2.gt.1.)
    call random_number(ran)
    xd=2.*(ran(1)-0.5)
    yd=2.*(ran(2)-0.5)
    zd=2.*(ran(3)-0.5)
    r2=xd**2+yd**2+zd**2
    !end do
    w=1./sqrt(r2)
    xd=xd*w
    yd=yd*w
    zd=zd*w
    !
  end subroutine random_versor
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !   place photon packets at a "star" surface
  !--------------------------------------------------------------------
  subroutine starsource(srad,is,js,ks,i,j,k,xd,yd,zd)
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
    rs= xs**2 + ys**2 + zs**2   
    w=float(srad)/sqrt(rs)
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
    rd = xd**2 + yd**2 + zd**2
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
  !   launches a photon from cell (xc,yc,zc) in the (xd,yd,zd) 
  !   direction, with an ionizing rate f
  !--------------------------------------------------------------------
  subroutine photons(xl0,yl0,zl0,xd,yd,zd,f)
    implicit none
    real, intent(in) :: xd,yd,zd
    real, intent(in) :: xl0, yl0, zl0
    real, intent(inout) :: f
    real :: xa, xb, ya, yb, za, zb
    real :: dl, dxl, dyl, dzl, xl, yl, zl, al, dtau
    real :: fmin
    integer :: iter, i, j ,k 
    !
    xa=1.
    xb=float(nx+1)
    ya=1.
    yb=float(ny+1)
    za=1.
    zb=float(nz+1)
    !
    !     photon trajectory
    !
    fmin=1.e-10*f
    dl=0.5
    dxl=dl*xd
    dyl=dl*yd
    dzl=dl*zd
    !
    xl=xl0+dxl
    yl=yl0+dyl
    zl=zl0+dzl
    !
    i=int(xl+0.5)
    j=int(yl+0.5)
    k=int(zl+0.5)  
    !
    
    !  the photon will cross the domain boundaries
    if (i < 1) then
       buffersize(1)=buffersize(1)+1
       photL(1,1, buffersize(1))= xl0+float(nx)
       photL(1,2, buffersize(1))= yl0
       photL(1,3, buffersize(1))= zl0
       !
       photL(1,4, buffersize(1))= xd
       photL(1,5, buffersize(1))= yd
       photL(1,6, buffersize(1))= zd
       photL(1,7, buffersize(1))= f
       return
    else if(i > nx) then
       buffersize(2)=buffersize(2)+1
       photR(1,1, buffersize(2))= xl0-float(nx)
       photR(1,2, buffersize(2))= yl0
       photR(1,3, buffersize(2))= zl0
       !
       photR(1,4, buffersize(2))= xd
       photR(1,5, buffersize(2))= yd
       photR(1,6, buffersize(2))= zd
       photR(1,7, buffersize(2))= f
       return
    endif
    if (j < 1) then
       buffersize(3)=buffersize(3)+1
       photB(1,1, buffersize(3))= xl0
       photB(1,2, buffersize(3))= yl0+float(ny)
       photB(1,3, buffersize(3))= zl0
       !
       photB(1,4, buffersize(3))= xd
       photB(1,5, buffersize(3))= yd
       photB(1,6, buffersize(3))= zd
       photB(1,7, buffersize(3))= f
       return
    else if(j > ny) then
       buffersize(4)=buffersize(4)+1
       photT(1,1, buffersize(4))= xl0
       photT(1,2, buffersize(4))= yl0-float(ny)
       photT(1,3, buffersize(4))= zl0
       !
       photT(1,4, buffersize(4))= xd
       photT(1,5, buffersize(4))= yd
       photT(1,6, buffersize(4))= zd
       photT(1,7, buffersize(4))= f
       return
    endif
    if (k < 1) then
       buffersize(5)=buffersize(5)+1
       photO(1,1, buffersize(5))= xl0
       photO(1,2, buffersize(5))= xl0
       photO(1,3, buffersize(5))= zl0+float(nz)
       !
       photO(1,4, buffersize(5))= xd
       photO(1,5, buffersize(5))= yd
       photO(1,6, buffersize(5))= zd
       photO(1,7, buffersize(5))= f
       return
    else if(k > nz) then
       buffersize(6)=buffersize(6)+1
       photI(1,1, buffersize(6))= xl0
       photI(1,2, buffersize(6))= yl0
       photI(1,3, buffersize(6))= zl0-float(nz)
       !
       photI(1,4, buffersize(6))= xd
       photI(1,5, buffersize(6))= yd
       photI(1,6, buffersize(6))= zd
       photI(1,7, buffersize(6))= f
       return
    end if
    !
    !
    do while( f.ge.fmin  &
         .and. (i <= nx).and.(i >= 1) &
         .and. (j <= ny).and.(k >= 1) &
         .and. (k <= nz).and.(k >= 1) &
         )
       !
       dtau=primit(neqdyn+1,i,j,k)*a0*dl*dx*rsc
       ph(i,j,k)=ph(i,j,k)+f*(1.-exp(-dtau) )/primit(neqdyn+1,i,j,k)
       f=f*exp(-dtau)
       !
       xl=xl+dxl
       yl=yl+dyl
       zl=zl+dzl
       !
       i=int(xl+0.5)
       j=int(yl+0.5)
       k=int(zl+0.5)
       !
       !  the photon reaches the domain boundaries
    if (i < 1) then
       buffersize(1)=buffersize(1)+1
       photL(1,1, buffersize(1))= xl-dxl+float(nx)
       photL(1,2, buffersize(1))= yl-dyl
       photL(1,3, buffersize(1))= zl-dzl
       !
       photL(1,4, buffersize(1))= xd
       photL(1,5, buffersize(1))= yd
       photL(1,6, buffersize(1))= zd
       photL(1,7, buffersize(1))= f
       return
    else if(i > nx) then
       buffersize(2)=buffersize(2)+1
       photR(1,1, buffersize(2))= xl-dxl-float(nx)
       photR(1,2, buffersize(2))= yl-dyl
       photR(1,3, buffersize(2))= zl-dzl
       !
       photR(1,4, buffersize(2))= xd
       photR(1,5, buffersize(2))= yd
       photR(1,6, buffersize(2))= zd
       photR(1,7, buffersize(2))= f
       return
    endif
    if (j < 1) then
       buffersize(3)=buffersize(3)+1
       photB(1,1, buffersize(3))= xl-dxl
       photB(1,2, buffersize(3))= yl-dyl+float(ny)
       photB(1,3, buffersize(3))= zl-dzl
       !
       photB(1,4, buffersize(3))= xd
       photB(1,5, buffersize(3))= yd
       photB(1,6, buffersize(3))= zd
       photB(1,7, buffersize(3))= f
       return
    else if(j > ny) then
       buffersize(4)=buffersize(4)+1
       photT(1,1, buffersize(4))= xl-dxl
       photT(1,2, buffersize(4))= yl-dyl-float(ny)
       photT(1,3, buffersize(4))= zl-dzl
       !
       photT(1,4, buffersize(4))= xd
       photT(1,5, buffersize(4))= yd
       photT(1,6, buffersize(4))= zd
       photT(1,7, buffersize(4))= f
       return
    endif
    if (k < 1) then
       buffersize(5)=buffersize(5)+1
       photO(1,1, buffersize(5))= xl-dxl
       photO(1,2, buffersize(5))= yl-dyl
       photO(1,3, buffersize(5))= zl-dzl+float(nz)
       !
       photO(1,4, buffersize(5))= xd
       photO(1,5, buffersize(5))= yd
       photO(1,6, buffersize(5))= zd
       photO(1,7, buffersize(5))= f
       return
    else if(k > nz) then
       buffersize(6)=buffersize(6)+1
       photI(1,1, buffersize(6))= xl-dxl
       photI(1,2, buffersize(6))= yl-dyl
       photI(1,3, buffersize(6))= zl-dzl-float(nz)
       !
       photI(1,4, buffersize(6))= xd
       photI(1,5, buffersize(6))= yd
       photI(1,6, buffersize(6))= zd
       photI(1,7, buffersize(6))= f
       return
    end if
       !
    end do
    !
    return
  end subroutine photons
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !  follows the rays across MPI boundaries    
  !--------------------------------------------------------------------
  subroutine radbounds()
#ifdef MPIP
    implicit none
    integer :: ip,jp,kp
    integer :: AllBufferSize(np*6), list(np), sizeSend, sizeRecv, niter
    integer,  parameter :: length=np*6
    integer:: status(MPI_STATUS_SIZE), err

    !   loop over MPI blocks to ensure the rays are followed to the 
    !   entire domain
    do ip=1,int( sqrt(float(mpicol)**2+float(mpirow)**2+float(mpirowz)**2) )

       !print'(a,2i3, 6i7)','**',rank,np,buffersize
       call mpi_allgather(buffersize(:)   , 6, mpi_integer, &
            Allbuffersize(:), 6, mpi_integer, comm2d ,err)

       !------------------------------------
       ! to the left
       if (left /= -1) then
          sizeSend=Allbuffersize(6*rank+1)
          if(sizeSend > 0) then
             call mpi_send(photL(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind    &
                  ,left,100,comm2d,err)
             !print*,rank,'send -> L:',SizeSend,'to   ',left
          end if
          !
          sizeRecv=Allbuffersize(6*left+2)
          if(sizeRecv > 0) then
             call mpi_recv(photL(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind   &
                  ,left,102, comm2d, status,err)
             !print*,rank,'recv <- L:',SizeRecv,'from ',left
          end if
       end if
       !------------------------------------
       !  to the right, step it up, step it up is allright
       if (right /= -1) then
          sizeRecv=Allbuffersize(6*right+1)
          if(sizeRecv > 0 ) then
             call mpi_recv(photR(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind  &
                  ,right,100,comm2d, status,err)
             !print*,rank, 'recv <- R:',SizeRecv,'from ',right
          end if
          sizeSend=Allbuffersize(6*rank +2)
          if(sizeSend > 0) then
             call mpi_send(photR(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind   &
                  ,right,102,comm2d,err)
             !print*,rank,'send -> R:',SizeSend,'to   ',right
          end if
       end if
       !
       !------------------------------------
       ! to the bottom
       if (bottom /= -1) then
          sizeSend=Allbuffersize(6*rank+3)
          if(sizeSend > 0) then
             call mpi_send(photB(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind    &
                  ,bottom,200,comm2d,err)
             !print*,rank,'send -> B:',SizeSend,'to   ',bottom
          end if
          !
          sizeRecv=Allbuffersize(6*bottom+4)
          if(sizeRecv > 0) then
             ! recv positions
             call mpi_recv(photB(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind   &
                  ,bottom,202, comm2d, status,err)
             !print*,rank,'recv <- B:',SizeRecv,'from ',bottom
          end if
       end if
       !
       !------------------------------------
       !  to the top
       if (top /= -1) then
          sizeRecv=Allbuffersize(6*top+3)
          if(sizeRecv > 0 ) then
             call mpi_recv(photT(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind  &
                  ,top,200,comm2d, status,err)
             print*,rank, 'recv <- T:',SizeRecv,'from ',top
          end if
          sizeSend=Allbuffersize(6*rank +4)
          if(sizeSend > 0) then
             call mpi_send(photT(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind   &
                  ,top,202,comm2d,err)
             !print*,rank,'send -> T:',SizeSend,'to   ',top
          end if
       end if
       !------------------------------------
       ! to the out
       if (out /= -1) then
          sizeSend=Allbuffersize(6*rank+5)
          if(sizeSend > 0) then
             call mpi_send(photO(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind    &
                  ,out,300,comm2d,err)
             !print*,rank,'send -> O:',SizeSend,'to   ',out
          end if
          !
          sizeRecv=Allbuffersize(6*out+6)
          if(sizeRecv > 0) then
             call mpi_recv(photO(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind   &
                  ,out,302, comm2d, status,err)
             !print*,rank,'recv <- O:',SizeRecv,'from ',out
          end if
       end if
       !------------------------------------
       !  to the in
       if (in /= -1) then
          sizeRecv=Allbuffersize(6*in+5)
          if(sizeRecv > 0 ) then
             call mpi_recv(photI(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind  &
                  ,in,300,comm2d, status,err)
             !print*,rank, 'recv <- I:',SizeRecv,'from ',in
          end if
          sizeSend=Allbuffersize(6*rank +6)
          if(sizeSend > 0) then
             call mpi_send(photI(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind   &
                  ,in,302,comm2d,err)
             !print*,rank,'send -> I:',SizeSend,'to   ',in
          end if
       end if
       !------------------------------------

       !  Continue with the photon tracing
       !  reset the out buffers
       buffersize(:)=0

       !  left face
       do niter=1,Allbuffersize(6*left+2)
          !print*,rank,'injecting photon',photL(2,:,niter)
          call photons(photL(2,1,niter),photL(2,2,niter),photL(2,3,niter), &
               photL(2,4,niter),photL(2,5,niter),photL(2,6,niter),photL(2,7,niter) )
       end do

       !  right face
       do niter=1,Allbuffersize(6*right+1)
          call photons(photR(2,1,niter),photR(2,2,niter),photR(2,3,niter), &
               photR(2,4,niter),photR(2,5,niter),photR(2,6,niter),photR(2,7,niter) )
       end do

       !  bottom face
       do niter=1,Allbuffersize(6*bottom+4)
          call photons(photB(2,1,niter),photB(2,2,niter),photB(2,3,niter), &
               photB(2,4,niter),photB(2,5,niter),photB(2,6,niter),photB(2,7,niter) )
       end do

       !  top face
       do niter=1,Allbuffersize(6*top+3)
          call photons(photT(2,1,niter),photT(2,2,niter),photT(2,3,niter), &
               photT(2,4,niter),photT(2,5,niter),photT(2,6,niter),photT(2,7,niter) )
       end do

       !  out face
       do niter=1,Allbuffersize(6*out+6)
          call photons(photO(2,1,niter),photO(2,2,niter),photO(2,3,niter), &
               photO(2,4,niter),photO(2,5,niter),photO(2,6,niter),photO(2,7,niter) )
       end do

       !  in face
       do niter=1,Allbuffersize(6*in+5)
          call photons(photI(2,1,niter),photI(2,2,niter),photI(2,3,niter), &
               photI(2,4,niter),photI(2,5,niter),photI(2,6,niter),photI(2,7,niter) )
       end do

       allBuffersize(:)=0

       !call mpi_barrier(mpi_comm_world, err)
       !call mpi_finalize(err)
       !stop

    end do
#endif
    return
  end subroutine radbounds
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
    real :: emax, f, dirx,diry,dirz, radius
    integer :: i,j,k, err, niter, ic,jc,kc, Srad, ii, jj, kk
    integer :: in, nmax

    !nmax = nxtot*nytot*nztot/10000/np
    !resets the photoionization rate
    ph (:,:,:)=0.
    em (:,:,:)=0.
    buffersize(:)=0
    !    
    !   computes the emissivity at each cell
    !call emdiff(emax)

    !   fire the photon torpedoes (nrays=1000000  moved to header)
    !   posicion de la fuente ionizante, in the entire domain
    ic=nxtot/2
    jc=nytot/2
    kc=nztot/2

    !   impose the photoionizing field of a star
    f=1.173e12*a0*1e3    !  flujo de fotones ionizantes:not divided by np
    srad=7

    in=0  ! number or rays successfully injected
    do niter=1, nrays 
       !  get the location and direction of the photon to be traced
       !  from 1:nxtot, 1:nytot, 1:nztot
       call starsource(srad,ic,jc,kc,i,j,k, dirx, diry, dirz)
       !  obtain in which proc will be traced
       ! hack
       ii=i/nx
       jj=j/ny
       kk=k/nz
       if( (ii==coords(0)).and.(jj==coords(1)).and.(kk==coords(2))) then
          in=in+1
          ! trace the photon
          call photons(float(i-coords(0)*nx),float(j-coords(1)*ny),float(k-coords(2)*nz), &
               dirx,diry,dirz,f)
          !call progress(niter,nrays)      
       end if

    end do

    !determine the actual number of photons injected, !
    !and divide ph among them
#ifdef MPIP
    call mpi_allreduce(in, nmax, 1, mpi_integer, mpi_sum, mpi_comm_world,err)
#else
    nmax=in
#endif
    ! trace photons across boundaries)
    call radbounds()
    !
    ph(:,:,:)=ph(:,:,:)/float(nmax)


    return
  end subroutine diffuse_rad
  !--------------------------------------------------------------------
end module difrad
#endif
