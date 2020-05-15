#
#  Makefile for Riemann solvers:
#
#  Type "make" or "make rp" to create .o files needed in 
#  applications directories
#
#
FFLAGS =  -c
LFLAGS =  
LINK    = gfortran
.f.o: ; $(LINK) $(FFLAGS) $*.f

RPOBJECTS = \
  rp1eu_HLLE.o 

RPSOURCES = \
  rp1eu_HLLE.f

rp: $(RPOBJECTS)


### DO NOT remove this line - make depends on it ###
