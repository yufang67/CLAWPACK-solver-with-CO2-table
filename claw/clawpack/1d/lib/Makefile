#
#		Makefile for the clawpack code:
#
#	To compile the claw1 library files type:  make lib
#       (from the AMRCLAWLIB sub-directory)
#
#	To make an executable, type:              make xclaw
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawhdf
#       (from the application sub-directory)
#
#	To compile a single file.f type:	  make file.o
#
#	To combine all source files type:	  make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
#	To combine all library source:    make backup
#       (date is appended)
#
FFLAGS = -O
LFLAGS =  
F77    = gfortran -c
LINK   = gfortran

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

LIBOBJECTS = \
  $(CLAW)/clawpack/1d/lib/claw1ez.o \
  $(CLAW)/clawpack/1d/lib/bc1.o \
  $(CLAW)/clawpack/1d/lib/setprob.o \
  $(CLAW)/clawpack/1d/lib/setaux.o \
  $(CLAW)/clawpack/1d/lib/b4step1.o \
  $(CLAW)/clawpack/1d/lib/claw1.o \
  $(CLAW)/clawpack/1d/lib/step1.o \
  $(CLAW)/clawpack/1d/lib/step1fw.o \
  $(CLAW)/clawpack/1d/lib/copyq1.o \
  $(CLAW)/clawpack/1d/lib/inlinelimiter.o \
  $(CLAW)/clawpack/1d/lib/limiter.o \
  $(CLAW)/clawpack/1d/lib/philim.o \
  $(CLAW)/clawpack/1d/lib/src1.o \
  $(CLAW)/clawpack/1d/lib/rp1.o \
  $(CLAW)/clawpack/1d/lib/out1.o \
  $(CLAW)/clawpack/1d/lib/out1_hdf.o

SOURCES = \
  $(CLAW)/clawpack/1d/lib/claw1ez.f \
  $(CLAW)/clawpack/1d/lib/bc1.f \
  $(CLAW)/clawpack/1d/lib/setprob.f \
  $(CLAW)/clawpack/1d/lib/setaux.f \
  $(CLAW)/clawpack/1d/lib/b4step1.f \
  $(CLAW)/clawpack/1d/lib/claw1.f \
  $(CLAW)/clawpack/1d/lib/step1.f \
  $(CLAW)/clawpack/1d/lib/step1fw.f \
  $(CLAW)/clawpack/1d/lib/copyq1.f \
  $(CLAW)/clawpack/1d/lib/inlinelimiter.f \
  $(CLAW)/clawpack/1d/lib/limiter.f \
  $(CLAW)/clawpack/1d/lib/philim.f \
  $(CLAW)/clawpack/1d/lib/src1.f \
  $(CLAW)/clawpack/1d/lib/rp1.f \
  $(CLAW)/clawpack/1d/lib/out1.f \
  $(CLAW)/clawpack/1d/lib/out1_hdf.f 

lib: $(LIBOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack1lib.`date +%y%m%d` 

clean:
	-rm -f $(LIBOBJECTS)

### DO NOT remove this line - make depends on it ###
