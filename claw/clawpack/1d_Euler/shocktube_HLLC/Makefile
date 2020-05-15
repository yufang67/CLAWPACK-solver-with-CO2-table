#
#               Makefile for the clawpack code:
#
#       To make an executable, type:            make xclaw
#       (from the application sub-directory)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To compile the lib library files type:    make lib
#       (from the claw/clawpack/1d/lib sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)

#
#
FFLAGS = -O
LFLAGS = -Wall -Wextra  -llapack -lblas
FCFLAGS += -I/$(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2 
F77    = gfortran -c
LINK   = gfortran

.f.o: ; $(F77) $(FFLAGS) $(FCFLAGS) $*.f -o $*.o

OBJECTS = \
  driver.o \
  qinit.o \
  rp_HLLC/rp1eu_HLLC.o \
  setprob.o 

LIBOBJECTS = \
  $(CLAW)/clawpack/1d/lib/setaux.o \
  $(CLAW)/clawpack/1d/lib/claw1ez.o \
  $(CLAW)/clawpack/1d/lib/bc1.o \
  $(CLAW)/clawpack/1d/lib/b4step1.o \
  $(CLAW)/clawpack/1d/lib/out1.o \
  $(CLAW)/clawpack/1d/lib/claw1.o \
  $(CLAW)/clawpack/1d/lib/step1.o \
  $(CLAW)/clawpack/1d/lib/copyq1.o \
  $(CLAW)/clawpack/1d/lib/limiter.o \
  $(CLAW)/clawpack/1d/lib/philim.o \
  $(CLAW)/clawpack/1d/lib/src1.o 

LIBTABLE = \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/def_constants.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/def_variables.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/deriv_disfonc.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/deriv_expfonc.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/helmholtz_deriv.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/helmholtz_dimless.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/properties.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/non_linear_solvers.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_functions.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_LL.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_LH.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_R.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_HT.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/saturation_curve.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/interp_functions.o \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/interp_table.o


LIBTABLESOURCES = \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/def_constants.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/def_variables.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/deriv_disfonc.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/deriv_expfonc.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/helmholtz_deriv.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/helmholtz_dimless.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/properties.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/non_linear_solvers.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_functions.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_LL.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_LH.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_R.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid_construction_HT.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/saturation_curve.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/grid.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/interp_functions.f90 \
  $(CLAW)/clawpack/1d_Euler/shocktube_HLLC/CO2/interp_table.f90


SOURCES = \
  driver.f \
  qinit.f \
  ../rp_HLLC/rp1eu_HLLC.f \
  setprob.f \
  $(CLAW)/clawpack/1d/lib/setaux.f \
  $(CLAW)/clawpack/1d/lib/claw1ez.f \
  $(CLAW)/clawpack/1d/lib/bc1.f \
  $(CLAW)/clawpack/1d/lib/b4step1.f \
  $(CLAW)/clawpack/1d/lib/out1.f \
  $(CLAW)/clawpack/1d/lib/claw1.f \
  $(CLAW)/clawpack/1d/lib/step1.f \
  $(CLAW)/clawpack/1d/lib/copyq1.f \
  $(CLAW)/clawpack/1d/lib/limiter.f \
  $(CLAW)/clawpack/1d/lib/philim.f \
  $(CLAW)/clawpack/1d/lib/src1.f 

xclaw: $(OBJECTS)
	$(LINK)  $(OBJECTS) $(LIBOBJECTS) $(LIBTABLE) -o xclaw $(LFLAGS)

program : $(SOURCES)
	cat $(SOURCES) > claw1program.f

clean:
	rm -f *.o *.mod *.MOD


### DO NOT remove this line - make depends on it ###
