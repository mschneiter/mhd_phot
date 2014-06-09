#---------------------------------------------------
#   for Huacho 3D
#---------------------------------------------------
#   Name of the executable
PROGRAM=PMHD_huacho
#
#   Choice of compiler and flags  ! if mpi mpif90 will be 
#   used but this is required to adjust some compilation flags
#COMPILER= gfortran
COMPILER= mpif90
#COMPILER= ifort
#FLAGS= -O3 -cpp -i_dynamic -mcmodel=large
#FLAGS= -O3 -cpp  -mcmodel=large -vec-report0 -traceback -check
#FLAGS= -O3 -cpp -vec-report0
FLAGS= -O3 -cpp
#
#---------------------------------------------------
#   compilation time parameters (Y=on, N=off)
#   carefull all of the following is case sensitive
#---------------------------------------------------
#   MPI paralalelization
MPIP = Y
#   Double Precision (Y/N Default is single)
DOUBLEP= Y
#
#   Enable Passive MHD 
#   Includes induction eq with no back reaction to the flow)
#   This is compatible with HD solvers (i.e. HLL, HLLC)   
PMHD = N
#
#   Enable MHD (UNDER CONSTRUCTION)
#   This is compatible with MHD solvers (HLLE, HLLD, under construction)
MHD = Y
#
#   Solver: 
#   HD  Solvers: HLL  (too difussive), or HLLC
#   MHD Solvers: HLLE (too difussive), or HLLD (under construction)
SOLVER = HLLC
#
#   Type of output
OUTDAT = N
OUTBIN = Y
OUTVTK = N
OUTSILO= N
#
#   IF silo was selected make sure the libraries are in place,
#   adjust the following line for that purpose
#
#FLAGS += -I/usr/local/silo/include -I/usr/local/szip/include -I/usr/local/hdf5/include/
#LINKFLAGS = -L/usr/local/szip/lib -L/usr/local/hdf5/lib/ -L/usr/local/silo/lib/
#LINKFLAGS += -lsiloh5  -lhdf5_fortran -lsz -lz -lstdc++ 
##   additional equations (i.e. passive scalars)?
PASSIVES = Y
#
#   Adiabatic or type of cooling (choose only one)
#   options are ADIABATIC, H, BBC, CORONAL
#   ADIABATIC: of course turns off the cooling
#   H: single parametrized cooling function (ionization fraction and T)
#   BBC: cooling function of Benjamin, Benson and Cox (2003)
#   DMC: Coronal eq. (tabulated) from Dalgarno & Mc Cray (1972)
#   CHI: from table generated with Chianti
COOLING = H
#
#   boundary conditions
#  (choose only one per boundary + OTHERB if needed)
#  optiins are OUTFLOW, CLOSED, PERIODIC
LEFTX   = OUTFLOW
RIGHTX  = OUTFLOW
BOTTOMY = OUTFLOW
TOPY    = OUTFLOW
INZ     = OUTFLOW
OUTZ    = OUTFLOW
OTHERB = N
#
#   choice of slope limiter, available limiters are:
#   limiter =-1: no average, 0 : no limiter, 1: Minmod,
#   2: Van Leer, 3: Van Albada,4: UMIST, 5: Woodward
#   6: Superbee
LIMITER =  1
#
#   Enable (isotropic) thermal conduction (Y/N)
THERMAL_COND = N
#
#   Enable diffuse radiation (Y/N) 
RADDIFF = Y
#  This version does the transfer locally
#  and passes the photons across mpi boundaries
#
#  Global diffuse radiation, is kept for testing purposes,
#  all the cores receive a copy of of the entire domain, does 
#  not allow very high resolution (deprecated)
RADDIFF_GLOBAL= N
#
# include source terms (Y/N)
SOURCE = Y
#
#   Include gravity (Y/N) (from point sources)
GRAV = N
#
#   Include radiative pressure (y/N)
RADPRES = N
#
#   Enable the N-Body module (Y/N) not implemented yet!
#   NBODY = N
#
#####################################################
# There should be no need to modify below this line #
#####################################################
#   List of objects
#---------------------------------------------------
OBJECTS = \
./src/parameters.o 	\
./src/globals.o	\
./src/dmc_module.o	\
./src/chi_module.o	\
./src/difrad.o \
./src/difrad_global.o \
./src/brio_wu.o \
./src/star.o \
./src/sources.o \
./src/thermal_cond.o \
./src/main.o		\
./src/initmain.o	\
./src/initflow.o 	\
./src/uprim.o 	\
./src/calcprim.o 	\
./src/Out_Silo_Module.o \
./src/output.o	\
./src/sound.o 	\
./src/timestep.o 	\
./src/tstep.o	\
./src/hllfluxes.o 	\
./src/hllcfluxes.o 	\
./src/hllefluxes.o 	\
./src/boundaryI.o	\
./src/boundaryII.o	\
./src/step.o 	\
./src/viscosity.o	\
./src/primf.o	\
./src/primu.o	\
./src/swaps.o	\
./src/limiter.o	\
./src/coolingh.o	\
./src/coolingdmc.o	\
./src/coolingchi.o	\
./src/bbcrd.o	\
./src/bbcev.o	\
./src/coolingbbc.o
#---------------------------------------------------
# ---- For the Lyman Alpha Tau calculation
OBJECTSLYAT = \
./src/parameters.o 	\
./src/globals.o	\
./src/uprim.o \
./src/lyman_alpha_tau.o \
#---------------------------------------------------
ifeq ($(DOUBLEP),Y) 
FLAGS += -DDOUBLEP
ifeq ($(COMPILER),ifort)
FLAGS += -r8
endif
ifeq ($(COMPILER),gfortran)
FLAGS += -fdefault-real-8
endif
endif
ifeq ($(SOLVER),HLL)
FLAGS += -DHLL
endif
ifeq ($(SOLVER),HLLC)
FLAGS += -DHLLC
endif
ifeq ($(SOLVER),HLLE)
FLAGS += -DHLLE
endif
ifeq ($(SOLVER),HLLD)
FLAGS += -DHLLD
endif
ifeq ($(OUTDAT),Y)
FLAGS += -DOUTDAT
endif
ifeq ($(OUTBIN),Y)
FLAGS += -DOUTBIN
endif
ifeq ($(OUTVTK),Y)
FLAGS += -DOUTVTK
endif
ifeq ($(OUTSILO),Y)
FLAGS += -DOUTSILO
endif
ifeq ($(PASSIVES),Y)
FLAGS += -DPASSIVES
endif
ifeq ($(COOLING),ADIABATIC)
FLAGS += -DADIABATIC
endif
ifeq ($(COOLING),H)
FLAGS += -DCOOLINGH
endif
ifeq ($(COOLING),DMC)
FLAGS += -DCOOLINGDMC
endif
ifeq ($(COOLING),CHI)
FLAGS += -DCOOLINGCHI
endif
ifeq ($(COOLING),BBC)
FLAGS += -DCOOLINGBBC
endif
ifeq ($(LEFTX),PERIODIC)
FLAGS += -DPERIODX
endif
ifeq ($(BOTTOMY),PERIODIC)
FLAGS += -DPERIODY
endif
ifeq ($(INZ),PERIODIC)
FLAGS += -DPERIODZ
endif
ifeq ($(LEFTX),CLOSED)
FLAGS += -DREFXL
endif
ifeq ($(RIGHTX),CLOSED)
FLAGS += -DREFXR
endif
ifeq ($(BOTTOMY),CLOSED)
FLAGS += -DREFYB
endif
ifeq ($(TOPY),CLOSED)
FLAGS += -DREFYT
endif
ifeq ($(INZ),CLOSED)
FLAGS += -DREFZI
endif
ifeq ($(OUTZ),CLOSED)
FLAGS += -DREFZO
endif
ifeq ($(LEFTX),OUTFLOW)
FLAGS += -DOUTFXL
endif
ifeq ($(RIGHTX),OUTFLOW)
FLAGS += -DOUTFXR
endif
ifeq ($(TOPY),OUTFLOW)
FLAGS += -DOUTFYT
endif
ifeq ($(BOTTOMY),OUTFLOW)
FLAGS += -DOUTFYB
endif
ifeq ($(INZ),OUTFLOW)
FLAGS += -DOUTFZI
endif
ifeq ($(OUTZ),OUTFLOW)
FLAGS += -DOUTFZO
endif
ifeq ($(OTHERB),Y)
FLAGS += -DOTHERB
endif
FLAGS += -DLIMITER=$(LIMITER)
ifeq ($(MPIP),Y)
FLAGS += -DMPIP
COMPILER = mpif90
endif
ifeq ($(RADDIFF),Y)
FLAGS += -DRADDIFF
endif
ifeq ($(RADDIFF_GLOBAL),Y)
FLAGS += -DRADDIFF_GLOBAL
endif
ifeq ($(NBODY),Y)
FLAGS += -DNBODY
endif
ifeq ($(THERMAL_COND),Y)
FLAGS += -DTHERMAL_COND
endif
ifeq ($(SOURCE),Y)
FLAGS += -DSOURCE
endif
ifeq ($(GRAV),Y)
FLAGS += -DGRAV
endif
ifeq ($(RADPRES),Y)
FLAGS += -DRADPRES
endif
ifeq ($(PMHD),Y)
FLAGS += -DPMHD
endif
#---------------------------------------------------
$(PROGRAM)  : ${OBJECTS}
	$(COMPILER) $(FLAGS) -o $@ $(OBJECTS) $(LINKFLAGS)

%.o:%.f95
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.f90
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.f
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.F95
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.F90
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.F
	$(COMPILER) $(FLAGS) -c $^ -o $@

.PHONY : clean
clean :
	rm -f *.o *.mod src/*.o src/*.mod
	rm -f $(PROGRAM) lyman_alpha_tau *.out *.err

cleanall :
	rm -f $(PROGRAM) lyman_alpha_tau *.out *.err
	rm -f *.o *.mod src/*.o src/*.mod
	rm -f BIN/*
	rm -f SILO/*.root
	rm -f SILO/BLOCKS/*

lyman_alpha_tau :

lyman_alpha_tau  : ${OBJECTSLYAT}
	$(COMPILER) $(FLAGS) -o $@ $(OBJECTSLYAT) $(LINKFLAGS)

%.o:%.f95
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.f90
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.f
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.F95
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.F90
	$(COMPILER) $(FLAGS) -c $^ -o $@
%.o:%.F
	$(COMPILER) $(FLAGS) -c $^ -o $@


#---------------------------------------------------
