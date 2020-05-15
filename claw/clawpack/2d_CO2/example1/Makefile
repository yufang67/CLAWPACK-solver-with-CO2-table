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
LFLAGS = $(FFLAGS)
F77    = f77 -c
LINK   = f77
HDFLIBS = -L/usr/local/hdf/lib -lmfhdf -ldf -ljpeg -lz


.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  driver.o \
  qinit.o \
  rpn2.o \
  rpt2.o \
  setprob.o 

LIBOBJECTS = \
  $(CLAW)/clawpack/2d/lib/claw2ez.o \
  $(CLAW)/clawpack/2d/lib/bc2.o \
  $(CLAW)/clawpack/2d/lib/setaux.o \
  $(CLAW)/clawpack/2d/lib/b4step2.o \
  $(CLAW)/clawpack/2d/lib/claw2.o \
  $(CLAW)/clawpack/2d/lib/step2.o \
  $(CLAW)/clawpack/2d/lib/step2ds.o \
  $(CLAW)/clawpack/2d/lib/dimsp2.o \
  $(CLAW)/clawpack/2d/lib/flux2.o \
  $(CLAW)/clawpack/2d/lib/copyq2.o \
  $(CLAW)/clawpack/2d/lib/limiter.o \
  $(CLAW)/clawpack/2d/lib/philim.o \
  $(CLAW)/clawpack/2d/lib/src2.o 

SOURCES = \
  driver.f \
  qinit.f \
  rpn2.f \
  rpt2.f \
  setprob.f \
  $(CLAW)/clawpack/2d/lib/claw2ez.f \
  $(CLAW)/clawpack/2d/lib/bc2.f \
  $(CLAW)/clawpack/2d/lib/setaux.f \
  $(CLAW)/clawpack/2d/lib/b4step2.f \
  $(CLAW)/clawpack/2d/lib/claw2.f \
  $(CLAW)/clawpack/2d/lib/step2.f \
  $(CLAW)/clawpack/2d/lib/step2ds.f \
  $(CLAW)/clawpack/2d/lib/dimsp2.f \
  $(CLAW)/clawpack/2d/lib/flux2.f \
  $(CLAW)/clawpack/2d/lib/copyq2.f \
  $(CLAW)/clawpack/2d/lib/limiter.f \
  $(CLAW)/clawpack/2d/lib/philim.f \
  $(CLAW)/clawpack/2d/lib/src2.f 

OUTOBJECTS = $(CLAW)/clawpack/2d/lib/out2.o \
	$(CLAW)/clawpack/2d/lib/restart2.o
OUTSOURCES = $(CLAW)/clawpack/2d/lib/out2.f \
	$(CLAW)/clawpack/2d/lib/restart2.f

HDFOBJECTS = $(CLAW)/clawpack/2d/lib/out2_hdf.o \
	$(CLAW)/clawpack/2d/lib/restart2_hdf.o
HDFSOURCES = $(CLAW)/clawpack/2d/lib/out2_hdf.f \
	$(CLAW)/clawpack/2d/lib/restart2_hdf.f

xclaw: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclaw

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
