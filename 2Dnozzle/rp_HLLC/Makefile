#
#  Makefile for Riemann solvers:
#
#  Type "make" or "make rp" to create .o files needed in 
#  applications directories
#
#
FFLAGS = -g -O0
LFLAGS = -c
FCFLAGS += -I/$(CLAW)/clawpack/CO2
F77    = gfortran -c
LINK   = gfortran

#.f.o: ; $(LINK)  $(LFLAGS) $(FCFLAGS) $*.f90 

%.o : %.f90
	$(LINK) $(LFLAGS) $(FCFLAGS)  $<


RPOBJECTS = \
  getquadinfo_mapgrid.o \
  euler_hllc_solver_mapgrid.o \
  rpt2_euler_mapgrid.o \
  rpn2_euler_mapgrid.o

RPSOURCES = \
  getquadinfo_mapgrid.f90 \
  euler_hllc_solver_mapgrid.f90 \
  rpt2_euler_mapgrid.f90 \
  rpn2_euler_mapgrid.f90

rp: $(RPOBJECTS) 

clean: 
	-rm -f  $(RPOBJECTS) 



### DO NOT remove this line - make depends on it ###
