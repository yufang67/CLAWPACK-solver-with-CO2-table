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
  $(CLAW)/clawpack/2d/lib/claw2ez.o \
  $(CLAW)/clawpack/2d/lib/bc2.o \
  $(CLAW)/clawpack/2d/lib/setprob.o \
  $(CLAW)/clawpack/2d/lib/setaux.o \
  $(CLAW)/clawpack/2d/lib/b4step2.o \
  $(CLAW)/clawpack/2d/lib/claw2.o \
  $(CLAW)/clawpack/2d/lib/step2.o \
  $(CLAW)/clawpack/2d/lib/step2ds.o \
  $(CLAW)/clawpack/2d/lib/dimsp2.o \
  $(CLAW)/clawpack/2d/lib/flux2.o \
  $(CLAW)/clawpack/2d/lib/flux2fw.o \
  $(CLAW)/clawpack/2d/lib/copyq2.o \
  $(CLAW)/clawpack/2d/lib/inlinelimiter.o \
  $(CLAW)/clawpack/2d/lib/limiter.o \
  $(CLAW)/clawpack/2d/lib/philim.o \
  $(CLAW)/clawpack/2d/lib/src2.o \
  $(CLAW)/clawpack/2d/lib/rpn2.o \
  $(CLAW)/clawpack/2d/lib/rpt2.o \
  $(CLAW)/clawpack/2d/lib/cellave.o \
  $(CLAW)/clawpack/2d/lib/fdisc.o \
  $(CLAW)/clawpack/2d/lib/out2.o \
  $(CLAW)/clawpack/2d/lib/out2_hdf.o \
  $(CLAW)/clawpack/2d/lib/restart2.o \
  $(CLAW)/clawpack/2d/lib/restart2_hdf.o

SOURCES = \
  $(CLAW)/clawpack/2d/lib/claw2ez.f \
  $(CLAW)/clawpack/2d/lib/bc2.f \
  $(CLAW)/clawpack/2d/lib/setprob.f \
  $(CLAW)/clawpack/2d/lib/setaux.f \
  $(CLAW)/clawpack/2d/lib/b4step2.f \
  $(CLAW)/clawpack/2d/lib/claw2.f \
  $(CLAW)/clawpack/2d/lib/step2.f \
  $(CLAW)/clawpack/2d/lib/step2ds.f \
  $(CLAW)/clawpack/2d/lib/dimsp2.f \
  $(CLAW)/clawpack/2d/lib/flux2.f \
  $(CLAW)/clawpack/2d/lib/flux2fw.f \
  $(CLAW)/clawpack/2d/lib/copyq2.f \
  $(CLAW)/clawpack/2d/lib/inlinelimiter.f \
  $(CLAW)/clawpack/2d/lib/limiter.f \
  $(CLAW)/clawpack/2d/lib/philim.f \
  $(CLAW)/clawpack/2d/lib/src2.f \
  $(CLAW)/clawpack/2d/lib/rpn2.f \
  $(CLAW)/clawpack/2d/lib/rpt2.f \
  $(CLAW)/clawpack/2d/lib/cellave.f \
  $(CLAW)/clawpack/2d/lib/fdisc.f \
  $(CLAW)/clawpack/2d/lib/out2.f \
  $(CLAW)/clawpack/2d/lib/out2_hdf.f \
  $(CLAW)/clawpack/2d/lib/restart2.f \
  $(CLAW)/clawpack/2d/lib/restart2_hdf.f

lib: $(LIBOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack2lib.`date +%y%m%d` 

clean:
	-rm -f $(LIBOBJECTS)

### DO NOT remove this line - make depends on it ###
