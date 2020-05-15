#
#		Makefile for the clawpack library code:
#
#	To compile the claw2 library files type:    make lib
#
#	To make an executable, type:                make xclaw
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:     make xclawhdf
#       (from the application sub-directory)
#
#	To compile a single file.f type:	    make file.o
#       (from the application sub-directory)
#
#	To combine source all source:               make backup
#       (date is appended)
#
FFLAGS = -O
LFLAGS =  
F77    = gfortran -c
LINK   = gfortran

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

LIBOBJECTS = \
  $(CLAW)/clawpack/2d_CO2/lib/claw2ez.o \
  $(CLAW)/clawpack/2d_CO2/lib/b4step2.o \
  $(CLAW)/clawpack/2d_CO2/lib/claw2.o \
  $(CLAW)/clawpack/2d_CO2/lib/step2ds.o \
  $(CLAW)/clawpack/2d_CO2/lib/dimsp2.o \
  $(CLAW)/clawpack/2d_CO2/lib/flux2fw.o \
  $(CLAW)/clawpack/2d_CO2/lib/flux2.o \
  $(CLAW)/clawpack/2d_CO2/lib/copyq2.o \
  $(CLAW)/clawpack/2d_CO2/lib/inlinelimiter.o \
  $(CLAW)/clawpack/2d_CO2/lib/limiter.o \
  $(CLAW)/clawpack/2d_CO2/lib/philim.o \
  $(CLAW)/clawpack/2d_CO2/lib/src2.o \
  $(CLAW)/clawpack/2d_CO2/lib/cellave.o \
  $(CLAW)/clawpack/2d_CO2/lib/fdisc.o \
  $(CLAW)/clawpack/2d_CO2/lib/out2_hdf.o \
  $(CLAW)/clawpack/2d_CO2/lib/restart2.o \
  $(CLAW)/clawpack/2d_CO2/lib/restart2_hdf.o

SOURCES = \
  $(CLAW)/clawpack/2d_CO2/lib/claw2ez.f \
  $(CLAW)/clawpack/2d_CO2/lib/b4step2.f \
  $(CLAW)/clawpack/2d_CO2/lib/claw2.f \
  $(CLAW)/clawpack/2d_CO2/lib/step2ds.f \
  $(CLAW)/clawpack/2d_CO2/lib/dimsp2.f \
  $(CLAW)/clawpack/2d_CO2/lib/flux2fw.f \
  $(CLAW)/clawpack/2d_CO2/lib/flux2.f \
  $(CLAW)/clawpack/2d_CO2/lib/copyq2.f \
  $(CLAW)/clawpack/2d_CO2/lib/inlinelimiter.f \
  $(CLAW)/clawpack/2d_CO2/lib/limiter.f \
  $(CLAW)/clawpack/2d_CO2/lib/philim.f \
  $(CLAW)/clawpack/2d_CO2/lib/src2.f \
  $(CLAW)/clawpack/2d_CO2/lib/cellave.f \
  $(CLAW)/clawpack/2d_CO2/lib/fdisc.f \
  $(CLAW)/clawpack/2d_CO2/lib/out2_hdf.f \
  $(CLAW)/clawpack/2d_CO2/lib/restart2.f \
  $(CLAW)/clawpack/2d_CO2/lib/restart2_hdf.f

lib: $(LIBOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack2lib.`date +%y%m%d` 

clean:
	-rm -f $(LIBOBJECTS)

### DO NOT remove this line - make depends on it ###
