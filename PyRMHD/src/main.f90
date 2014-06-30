!=======================================================================
!   integrates Euler equations in three dimensions
!   the choice of the integration method is set in the makefile
!
!   the flow variables:
!   ieq=1 : rho  (total)
!       2 : rho u
!       3 : rho v
!       4 : rho w
!       5 : Internal energy (thermal+kinetic)
!       6 : bx
!       7 : by
!       8 : bz
!
!       additional variables advected into the flow, e.g.:
!       9  (6):  n_HI
!       10 (7):  n_HII
!       8  (8):  n_HeI
!       9  (9):  n_HeII
!      10  (10): n_HeIII
!      11  (11): rho*zbar
!      12  (12): ne
!      But they can be changed bu the user according to cooling
!      function for instance
!=======================================================================
program euler
  use parameters
  use globals
#ifdef RADDIFF
  use difrad
#endif
#ifdef RADDIFF_GLOBAL
  use difrad_global
#endif
  implicit none
  integer :: err
  integer :: itprint
  real    :: time, dt, tprint
  !--------------------------------------------------------------------
  !   initializes mpi, and global variables
  call initmain(time, tprint, itprint)
  !---------------------------------------------------------------------
  !   set initial conditions
  !   initialize u's
  call initflow(itprint)

  !---------------------------------------------------------------------
  !   impose  boundaries (needed if imposing special BCs)
  call boundaryI(time,0.0)

  !---------------------------------------------------------------------
  !   update primitives with u
  call calcprim(u,primit)

  !---------------------------------------------------------------------
  !  if the radiative transfer is enabled, proceed to it
#ifdef RADDIFF
  call diffuse_rad()
#endif
#ifdef RADDIFF_GLOBAL
  call diffuse_rad()
#endif

#ifdef MPIP
  call mpi_barrier(mpi_comm_world, err)
#endif

  !------------------------------------------------------------------
  !   time integration
  do while (time.le.tmax)
     !-----------------------------------
     !   output at intervals tprint
     if(time.ge.tprint) then
        !call diffuse_rad()
        call output(itprint)
        if (rank.eq.0) then
           print'(a,i4)','****************** wrote output *************** &
                & :' , itprint
        end if
        tprint=tprint+dtprint
        itprint=itprint+1
     end if
     !-----------------------------------

     !   computes the timestep
     call timestep(dt)
     !if ((dt*tsc/yr) .gt. 1000.) then
     ! dt=1.E3*yr/tsc
     !endif
     time = time + dt
     if (rank.eq.1) print'(a,es12.3,a,es12.3,a,es12.3,a)',            &
     'time=',time,'  dt=', dt,' tprint=',tprint
     !'time=',time*tsc/day,' day  dt=', dt*tsc/day,' day tprint=',tprint*tsc/day,' day'
     !'time=',time*tsc/yr,' yr  dt=', dt*tsc/yr,' yr tprint=',tprint*tsc/yr,' yr'
     !
     !   advances the solution
     call tstep(time,dt)
     !
  end do
  !   finishes
  if (rank.eq.0) print'(a)',"--- My work here is done, have a nice day ---"
  !---------------------------------------------------------------------
#ifdef MPIP
  call mpi_finalize(err)
#endif
  stop
end program  euler
!=======================================================================
