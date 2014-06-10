!=======================================================================
!  initialization routine
!=======================================================================
subroutine initmain(time, tprint, itprint)
  use parameters
  use globals
#ifdef COOLINGDMC
  use dmc_module
#endif
#ifdef COOLINGCHI
  use chi_module
#endif
#ifdef RADDIFF
  use difrad
#endif
#ifdef RADDIFF_GLOBAL
  use difrad_global
#endif
#ifdef THERMAL_COND
  use thermal_cond
#endif
  use user_mod

  implicit none
  real,    intent(out) :: time, tprint
  integer, intent(out) :: itprint
#ifdef MPIP
  integer :: ilevs, err, nps
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period
#endif
#ifdef COOLINGDMC
  integer :: i
  real (kind=8) :: a, b
#endif
#ifdef COOLINGCHI
  integer :: i
  real (kind=8) :: a, b
#endif
  !
  !--------------------------------------------------------------------
  !initializes MPI
#ifdef MPIP
#ifdef PERIODX
  logical, parameter :: perx=.true.
#else
  logical, parameter :: perx=.false.
#endif
#ifdef PERIODY
  logical, parameter :: pery=.true.
#else
  logical, parameter :: pery=.false.
#endif
#ifdef PERIODZ
  logical, parameter :: perz=.true.
#else
  logical, parameter :: perz=.false.
#endif
  period(0)=perx
  period(1)=pery
  period(2)=perz
  dims(0)  =mpicol
  dims(1)  =mpirow
  dims(2)  =mpirowz
  !
  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)
  if (nps.ne.np) then
     print*, 'processor number (',nps,') is not equal to pre-defined number (',np,')'
     call mpi_finalize(err) 
     stop
  endif
#else
  rank=0
  coords(:)=0
#endif
  if(rank.eq.master) then
     print '(a)' ,"*******************************************"
     print '(a)' ," _                      _                 *"
     print '(a)' ,"| |__  _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ,"| '_ \| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| | | | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ,"|_| |_|\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ,"                                          *"
  endif
#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', np , ' processors    *'
     print '(a)' ,'*******************************************'
  end if
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, 1            &
       , comm2d, err)
  call mpi_comm_rank(comm2d, rank, err)
  call mpi_cart_coords(comm2d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                              &
       ,' ready w/coords',coords(0),coords(1),coords(2)   
  call mpi_cart_shift(comm2d, 0, 1, left  , right, err)
  call mpi_cart_shift(comm2d, 1, 1, bottom, top  , err)
  call mpi_cart_shift(comm2d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)   
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
#endif
  !
  !   grid spacing
  dx=xmax/nxtot
  dy=ymax/nytot
  dz=zmax/nztot
  !
  !   initialize time integration
  if (iwarm.eq.0) then
     if(rank.eq.master) then
        print'(a)', 'Starting cold'
        print'(a)',' ' 
     endif
     itprint=0
     time=0.
     tprint=0.
  else
     itprint=itprint0
     time=itprint*dtprint
!     if(rank.eq.master) then
!        print'(a,i,a,es12.3,a)', 'Warm start , from output ',itprint,' at a tim!e ',time*tsc/yr,' yr'
!        print'(a)',' ' 
!     end if
     tprint=time+dtprint
  end if
  !
  !   allocate big arrays in memory
  allocate (     u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (    up(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (primit(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (     f(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (     g(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (     h(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  !
  !   DMC cooling
#ifdef COOLINGDMC
  if(rank.eq.master) then
     open(unit=10,file=trim(workdir)//'/source/DMClib/coolingDMC.tab',status='old')
     do i=1,41
        read(10,*) a, b
        cooltab(1,i)=10.d0**(a)
        cooltab(2,i)=10.d0**(-b)
     end do
     close(unit=10)
  endif
#ifdef MPIP
  call mpi_bcast(cooltab,82,mpi_double_precision,0,mpi_comm_world,err)
#endif
#endif
!   CHIANTI COOLING
#ifdef COOLINGCHI
  if(rank.eq.master) then
     open(unit=10,file= trim(workdir)//'/source/CHIANTIlib/coolingCHIANTI.tab',status='old')
     do i=1,41
        read(10,*) a, b
        cooltab(1,i)=a
        cooltab(2,i)=b
     end do
     close(unit=10)
  endif
#ifdef MPIP
  call mpi_bcast(cooltab,82,mpi_double_precision,0,mpi_comm_world,err)
#endif
#endif
!   BBC COOLING
#ifdef COOLINGBBC
  do ii=0,(nps-1)
     if (rank.eq.ii) then
        call bbcrd
        print'(a,i4,a)','rank:',rank,' Just read the tables'
     endif
#ifdef MPIP
     call mpi_barrier (mpi_comm_world, err)
#endif
  end do
#endif
  !
#ifdef THERMAL_COND
  call init_thermal_cond()
#endif
#ifdef RADDIFF
  call init_rand()
#endif
#ifdef RADDIFF_GLOBAL
  call init_rand()
#endif
  !  User input initialization, it is called always, 
  !  it has to be there, even empty
  call init_user_mod()
  !
  !--------------------------------------------------------------------
  !   write report of compilation parameters
#ifdef MPIP
  call mpi_barrier(mpi_comm_world, err)
  if(rank.eq.master) then
     print'(a)',''
#endif
#ifdef DOUBLEP
     print'(a)', 'Double precision used (reals are 8 bytes long)'
     print'(a)', ''
#else
     print'(a)', 'Single precision used (reals are 4 bytes long)'
     print'(a)', ''
#endif
#ifdef HLLC
     print'(a)', 'The Riemann solver is HLLC'
     print'(a)', ''
#endif
#ifdef HLL
     print'(a)', 'The Riemann solver is HLL'
     print'(a)', ''
#endif
#ifdef HLL_HLLC
     print'(a)', 'The Riemann solver is HLL-HLLC (hybrid)'
     print'(a)', ''
#endif
#ifdef ADIABATIC
     print'(a)', 'The code is in ADIABATIC mode'
     print'(a)', ''
#endif
#ifdef COOLINGH
     print'(a)', 'Radiative cooling ON (w/parametrized cooling curve)'
     print'(a)', ''
#endif
#ifdef COOLINGBBG
     print'(a)', 'Radiative cooling ON (w/ Benjamin Benson & Cox 2003 prescription)'
     print'(a)', ''
#endif
#ifdef COOLINGDMC
     print'(a)', 'Radiative cooling ON (w/ Dalgarno & Mc Cray, coronal eq.)'
     print'(a)', ''
#endif
#ifdef COOLINGCHI
     print'(a)', 'Radiative cooling ON (Uses table from CHIANTI)'
     print'(a)', ''
#endif
#ifdef RADDIFF_GLOBAL
     print'(a)','Diffuse radiative transfer enabled, global'
#endif
#ifdef RADDIFF
     print'(a)','Diffuse radiative transfer enabled, local'
#endif
     print'(a)', '-----  OUTPUT -----------------------'
     print'(a)', 'path: '//trim(outputpath)
     print'(a)', 'in the following format(s):'
#ifdef OUTBIN
     print'(a)', '*.bin (raw unformatted)'
#endif
#ifdef OUTDAT
     print'(a)', '*.dat (formatted, beware of big files)'
#endif
#ifdef OUTVTK
     print'(a)', '*.vtk (binary VTK)'
#endif
     print'(a)', ''
     print'(a)', '----- BOUNDARY CONDITIONS -----------'
#ifdef PERIODX
     print'(a)', 'LEFT & RIGHT: PERIODIC'
#endif
#ifdef PERIODY
     print'(a)', 'BOTTOM & TOP: PERIODIC'
#endif
#ifdef PERIODZ
     print'(a)', 'IN & OUT: PERIODIC'
#endif
#ifdef REFXL
     print'(a)', 'LEFT:   REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFXL
     print'(a)', 'LEFT:   OUTFLOW    (OPEN)'
#endif
#ifdef INFXL
     print'(a)', 'LEFT:   INFLOW     (user defined)'
#endif
#ifdef REFXR
     print'(a)', 'RIGHT:  REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFXR
     print'(a)', 'RIGHT:  OUTFLOW    (OPEN)'
#endif
#ifdef INFXR
     print'(a)', 'RIGHT:  INFLOW     (user defined)'
#endif
#ifdef REFYB
     print'(a)', 'BOTTOM: REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFYB
     print'(a)', 'BOTTOM: OUTFLOW    (OPEN)'
#endif
#ifdef INFYB
     print'(a)', 'BOTTOM: INFLOW     (user defined)'
#endif
#ifdef REFYT
     print'(a)', 'TOP:    REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFYT
     print'(a)', 'TOP:    OUTFLOW    (OPEN)'
#endif
#ifdef INFYT
     print'(a)', 'TOP:    INFLOW     (user defined)'
#endif
#ifdef REFZO
     print'(a)', 'OUT:    REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFZO
     print'(a)', 'OUT:    OUTFLOW    (OPEN)'
#endif
#ifdef INFZO
     print'(a)', 'OUT:    INFLOW     (user defined)'
#endif
#ifdef REFZI
     print'(a)', 'IN: REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFZI
     print'(a)', 'IN:     OUTFLOW    (OPEN)'
#endif
#ifdef INFZI
     print'(a)', 'IN: INFLOW     (user defined)'
#endif
    print'(a)', ''
     print'(a)', '----- OTHER STUFF -----------'
#ifdef RADDIFF
     print'(a)', 'Diffuse radiation (local+MPI) enabled'
#endif
#ifdef RADDIFF_GLOBAL
     print'(a)', 'Diffuse radiation (global, RAM hungry) enabled'
#endif
#ifdef THERMAL_COND
     print'(a)', 'Thermal conduction included (isotropic)'
#endif
#ifdef OTHERB
     print'(a)', 'Other boundaries enabled (otherbounds.f90)'
#endif
     print'(a)', ''
#if LIMITER==-1
     print'(a)', 'No average in the limiter (reduces to 1st order)'
#endif
#if LIMITER==0
     print'(a)', 'No limiter'
#endif
#if LIMITER==1
     print'(a)', 'MINMOD limiter -most diffusive-'
#endif
#if LIMITER==2
     print'(a)', 'Falle Limiter (Van Leer)'
#endif
#if LIMITER==3
     print'(a)', 'Van Albada Limiter'
#endif
#if LIMITER==4
     print'(a)', 'UMIST limiter -least diffusive-'
#endif
#if LIMITER==5
     print'(a)', 'Woodward Limiter (MC-limiter; monotonized central difference)'
#endif
#if LIMITER==6
     print'(a)', 'SUPERBEE limiter (tends to flatten circular waves)'
#endif
     print'(a)', ''
     print'(a)','***********************************************'
#ifdef MPIP
  end if
  call mpi_barrier(mpi_comm_world, err)
#endif
  !--------------------------------------------------------------------
end subroutine initmain
  !====================================================================
