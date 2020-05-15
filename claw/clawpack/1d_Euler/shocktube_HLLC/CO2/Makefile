#
#  Makefile for CO2 table:
#
#  Type "make" or "make rp" to create .o files needed in 
#  applications directories
#
#
#
FFLAGS = -g -O0
LFLAGS = -c
F77    = gfortran -c
LINK   = gfortran

#.f.o: ; $(LINK)  $(LFLAGS) $*.f90 
 
%.o : %.f90
	$(LINK) $(LFLAGS)  $<

OBJECTS = \
  def_constants.o \
  def_variables.o \
  deriv_disfonc.o \
  deriv_expfonc.o \
  helmholtz_deriv.o \
  helmholtz_dimless.o \
  properties.o \
  non_linear_solvers.o \
  grid_functions.o \
  grid_construction_LL.o \
  grid_construction_LH.o \
  grid_construction_R.o \
  grid_construction_HT.o \
  saturation_curve.o \
  grid.o \
  interp_functions.o \
  interp_table.o
 

CO2SOURCES = \
  def_constants.f90 \
  def_variables.f90 \
  deriv_disfonc.f90 \
  deriv_expfonc.f90 \
  helmholtz_deriv.f90 \
  helmholtz_dimless.f90 \
  properties.f90 \
  non_linear_solvers.f90 \
  grid_functions.f90 \
  grid_construction_LL.f90 \
  grid_construction_LH.f90 \
  grid_construction_R.f90 \
  grid_construction_HT.f90 \
  saturation_curve.f90 \
  grid.f90 \
  interp_functions.f90 \
  interp_table.f90

CO2: $(OBJECTS)

clean:
	rm -f *.o *.mod *.MOD



### DO NOT remove this line - make depends on it ###
