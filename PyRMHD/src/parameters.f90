!=======================================================================
!   This module contains all variables and constants that are global
!=======================================================================
module parameters
  implicit none
#ifdef MPIP
  include "mpif.h"
#endif
  !
 character (len=128),parameter ::  outputpath='./' &
                                   ,  workdir='./'
  !
#if defined(MHD) || defined(PMHD)
  integer, parameter :: neqdyn=8        !# of eqs  (+scal)
#else
  integer, parameter :: neqdyn=5        !# of eqs  (+scal)
#endif
  integer, parameter :: ndim=3          !dimensions
  integer, parameter :: npas=0          ! number of passive scalars 
  integer, parameter :: nghost=2        ! number of ghost cells
  !   mesh size
  integer, parameter :: nxtot=400
  integer, parameter :: nytot=8
  integer, parameter :: nztot=8
  !
#ifdef MPIP
  !   mpi array of processors
  integer, parameter :: mpicol=8, mpirow=1,mpirowz=1
  integer, parameter :: np=mpicol*mpirow*mpirowz
#endif
  !---------------------------------------------------------------------
  !   the following options could go on an runtime input file
  !---------------------------------------------------------------------
  !
  !   some physical constants (and pi)
  real, parameter :: amh=1.66e-24,mu=1.3,Kb=1.38e-16
  real, parameter :: Hmm=1.007,Rg=8.3145e7, Ggrav=6.67259e-8
  real, parameter :: cv=1., gamma=(cv+1.)/cv
  real, parameter :: pi=3.1415926535
  real, parameter :: Msun=1.99E33, Rsun=6.955e10
  real, parameter :: Mjup=1.898E30
  real, parameter :: pc=3.0857E18, AU=1.496e13
  real, parameter :: clight=2.99E10
  !  some other constants
  real, parameter :: yr=3.1536E7, day=86400., hr=3600. !yr, day and hour in sec.
  !
  !   box size in code units
  real, parameter :: xmax=1., ymax=0.01, zmax=0.01
  !   box size, corresponds to xmax
  real, parameter :: xphys=1.
  !
  !---------------------------------------------------------------------
  !
  !   scaling factors to physical (cgs) units
  real, parameter :: T0=1.e4
  real, parameter :: rsc=xphys/xmax
  real, parameter :: rhosc=amh*mu 
  real, parameter :: Tempsc=T0*gamma   
  real, parameter :: vsc2 =gamma*Rg*T0/mu
  real, parameter :: Psc = rhosc*vsc2
  real, parameter :: tsc =rsc/sqrt(vsc2) 
#ifdef PMHD
   real, parameter :: bsc = sqrt(4.0*pi*Psc)     
#endif
#ifdef MHD
   real, parameter :: bsc = sqrt(4.0*pi*Psc)     
#endif
  !----------------------------------------------------
  !   flow properties of initial conditions (cgs)
  !   environment (cgs) / SN  parameters are in initconds
  real, parameter :: nenv=44., Tenv=10., venv=0.0
  !--------------------------------------------------------------------
  !   maximum time of integration and output interval
  real, parameter :: tmax    =0.15 !*day/tsc
  real, parameter :: dtprint =0.02 !*day/tsc
  real, parameter :: cfl=0.5                    ! Courant number
  real, parameter :: eta=0.001                  ! viscosity
  !
  !   for an iwarm (0->from t=0,1->starts from t=itprint*dtprint
  integer, parameter :: iwarm=0
  integer            :: itprint0=4
  !--------------------------------------------------------------------
  !   some derived parameters (no need of user's input below this line)
  !--------------------------------------------------------------------
  !
  integer, parameter :: neq=neqdyn + npas
!  integer, parameter :: neqpas= neq
#ifdef MPIP
  integer, parameter :: nx=nxtot/mpicol, ny=nytot/mpirow, nz=nztot/mpirowz
#else
  integer, parameter :: nx=nxtot, ny=nytot, nz=nztot
  integer, parameter :: np=1, mpirow=1, mpicol=1, mpirowz=1
#endif
  !  array bounds 
  integer, parameter :: nxmin=1-nghost, nxmax=nx+nghost
  integer, parameter :: nymin=1-nghost, nymax=ny+nghost
  integer, parameter :: nzmin=1-nghost, nzmax=nz+nghost 
  !
  !   more mpi stuff
  integer, parameter ::master=0
  !
  !   set floating point precision (kind) for MPI messages
#ifdef MPIP
#ifdef DOUBLEP
  integer, parameter :: mpi_real_kind=mpi_real8
#else
  integer, parameter :: mpi_real_kind=mpi_real4
#endif
#endif
  !---------------------------------------------------------------------
end module parameters
!=======================================================================

