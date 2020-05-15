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
.f.o: ; $(F77) $(FFLAGS) $(FCFLAGS) $*.f -o $*.o
#%.o : %.f90
#	$(LINK) $(LFLAGS) $(FCFLAGS)  $<


RPOBJECTS = \
  rp1eu_HLLC.o  

RPSOURCES = \
  rp1eu_HLLC.f 

rp: $(RPOBJECTS) 

clean: 
	-rm -f  $(RPOBJECTS) 



### DO NOT remove this line - make depends on it ###
