#
#               Makefile for the clawpack code:
#
#       To make an executable, type:            make xclaw
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawhdf
#       (from the application sub-directory)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To compile the lib library files type:    make lib
#       (from the claw/clawpack/2d/lib sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
FFLAGS = -O
LFLAGS = -Wall -Wextra  -llapack -lblas
FCFLAGS += -I/$(CLAW)/clawpack/CO2
F77    = gfortran -c
LINK   = gfortran
HDFLIBS = -L/usr/local/hdf/lib -lmfhdf -ldf -ljpeg -lz

.f.o: ; $(F77) $(FFLAGS) $(FCFLAGS) $*.f -o $*.o




OBJECTS = \
  driver.o \
  qinit.o \
  setprob.o \
  mapc2p.o \
  setaux.o \
  out2.o \
  flux2.o \
  step2.o \
  restart2.o \
  step2ds.o \
  dimsp2.o \
  bc2.o

LIBOBJECTS = \
  $(CLAW)/clawpack/2d_CO2/lib/claw2ez.o \
  $(CLAW)/clawpack/2d_CO2/lib/b4step2.o \
  $(CLAW)/clawpack/2d_CO2/lib/claw2.o \
  $(CLAW)/clawpack/2d_CO2/lib/copyq2.o \
  $(CLAW)/clawpack/2d_CO2/lib/limiter.o \
  $(CLAW)/clawpack/2d_CO2/lib/philim.o \
  $(CLAW)/clawpack/2d_CO2/lib/src2.o \
  rp_HLLC/euler_hllc_solver_mapgrid.o \
  rp_HLLC/rpt2_euler_mapgrid.o \
  rp_HLLC/rpn2_euler_mapgrid.o \
  rp_HLLC/getquadinfo_mapgrid.o


SOURCES = \
  driver.f \
  qinit.f \
  setprob.f \
  mapc2p.f\
  setaux.f\
  bc2.f \
  out2.f \
  flux2.f \
  step2.f \
  restart2.f \
  ./rp_HLLC/euler_hllc_solver_mapgrid.f90 \
  ./rp_HLLC/rpt2_euler_mapgrid.f90 \
  ./rp_HLLC/rpn2_euler_mapgrid.f90 \
  ./rp_HLLC/getquadinfo_mapgrid.f90 \
  $(CLAW)/clawpack/2d_CO2/lib/claw2ez.f \
  $(CLAW)/clawpack/2d_CO2/lib/b4step2.f \
  $(CLAW)/clawpack/2d_CO2/lib/claw2.f \
  step2ds.f \
  dimsp2.f \
  $(CLAW)/clawpack/2d_CO2/lib/copyq2.f \
  $(CLAW)/clawpack/2d_CO2/lib/limiter.f \
  $(CLAW)/clawpack/2d_CO2/lib/philim.f \
  $(CLAW)/clawpack/2d_CO2/lib/src2.f 


HDFOBJECTS = $(CLAW)/clawpack/2d_CO2/lib/out2_hdf.o \
	$(CLAW)/clawpack/2d_CO2/lib/restart2_hdf.o
HDFSOURCES = $(CLAW)/clawpack/2d_CO2/lib/out2_hdf.f \
	$(CLAW)/clawpack/2d_CO2/lib/restart2_hdf.f

LIBTABLE = \
  $(CLAW)/clawpack/CO2/Def_constants.o \
  $(CLAW)/clawpack/CO2/Def_variables.o \
  $(CLAW)/clawpack/CO2/deriv_disfonc.o \
  $(CLAW)/clawpack/CO2/deriv_expfonc.o \
  $(CLAW)/clawpack/CO2/helmholtz_deriv.o \
  $(CLAW)/clawpack/CO2/helmholtz_dimless.o \
  $(CLAW)/clawpack/CO2/Properties.o \
  $(CLAW)/clawpack/CO2/non_linear_solvers.o \
  $(CLAW)/clawpack/CO2/grid_functions.o \
  $(CLAW)/clawpack/CO2/grid_construction_LL.o \
  $(CLAW)/clawpack/CO2/grid_construction_LH.o \
  $(CLAW)/clawpack/CO2/grid_construction_R.o \
  $(CLAW)/clawpack/CO2/grid_construction_HT.o \
  $(CLAW)/clawpack/CO2/grid_construction_TPL.o \
  $(CLAW)/clawpack/CO2/grid_construction_TPM.o \
  $(CLAW)/clawpack/CO2/grid_construction_TPH.o \
  $(CLAW)/clawpack/CO2/saturation_curve.o \
  $(CLAW)/clawpack/CO2/Grid.o \
  $(CLAW)/clawpack/CO2/interp_functions.o \
  $(CLAW)/clawpack/CO2/Derivees.o \
  $(CLAW)/clawpack/CO2/Interp_table.o \
  $(CLAW)/clawpack/CO2/peng_robinson.o \
  $(CLAW)/clawpack/CO2/stiffened.o \
  $(CLAW)/clawpack/CO2/solver_eos.o \
  $(CLAW)/clawpack/CO2/var_const.o \
  $(CLAW)/clawpack/CO2/location.o


LIBTABLESOURCES = \
  $(CLAW)/clawpack/CO2/Def_constants.f90 \
  $(CLAW)/clawpack/CO2/Def_variables.f90 \
  $(CLAW)/clawpack/CO2/deriv_disfonc.f90 \
  $(CLAW)/clawpack/CO2/deriv_expfonc.f90 \
  $(CLAW)/clawpack/CO2/helmholtz_deriv.f90 \
  $(CLAW)/clawpack/CO2/helmholtz_dimless.f90 \
  $(CLAW)/clawpack/CO2/Properties.f90 \
  $(CLAW)/clawpack/CO2/non_linear_solvers.f90 \
  $(CLAW)/clawpack/CO2/grid_functions.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_LL.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_LH.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_R.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_HT.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_TPH.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_TPL.f90 \
  $(CLAW)/clawpack/CO2/grid_construction_TPM.f90 \
  $(CLAW)/clawpack/CO2/saturation_curve.f90 \
  $(CLAW)/clawpack/CO2/Grid.f90 \
  $(CLAW)/clawpack/CO2/interp_functions.f90 \
  $(CLAW)/clawpack/CO2/Derivees.f90 \
  $(CLAW)/clawpack/CO2/Interp_table.f90 \
  $(CLAW)/clawpack/CO2/peng_robinson.f90 \
  $(CLAW)/clawpack/CO2/stiffened.f90 \
  $(CLAW)/clawpack/CO2/var_const.f90 \
  $(CLAW)/clawpack/CO2/solver_eos.f90 
                                           

xclaw: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(LIBTABLE) -o xclaw

xclawhdf: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) $(HDFLIBS) \
	-o xclawhdf

program : $(SOURCES) $(OUTSOURCES)
	cat $(SOURCES) $(OUTSOURCES) > claw2program.f

programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > claw2program.f

clean:
	-rm -f $(OBJECTS) xclaw xclawhdf 

### DO NOT remove this line - make depends on it ###
