#
#               Makefile for the clawpack code:
#
#       To make an executable, type:              make xclaw
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
#       (from the claw/clawpack/3d/lib sub-directory)
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
  rpn3bu.o \
  rpt3bu.o \
  rptt3bu.o \
  setprob.o

LIBOBJECTS = \
  $(CLAW)/clawpack/3d/lib/claw3ez.o \
  $(CLAW)/clawpack/3d/lib/bc3.o \
  $(CLAW)/clawpack/3d/lib/setaux.o \
  $(CLAW)/clawpack/3d/lib/b4step3.o \
  $(CLAW)/clawpack/3d/lib/chkmth.o \
  $(CLAW)/clawpack/3d/lib/claw3.o \
  $(CLAW)/clawpack/3d/lib/step3.o \
  $(CLAW)/clawpack/3d/lib/step3ds.o \
  $(CLAW)/clawpack/3d/lib/dimsp3.o \
  $(CLAW)/clawpack/3d/lib/flux3.o \
  $(CLAW)/clawpack/3d/lib/copyq3.o \
  $(CLAW)/clawpack/3d/lib/limiter.o \
  $(CLAW)/clawpack/3d/lib/philim.o \
  $(CLAW)/clawpack/3d/lib/src3.o 

SOURCES = \
  driver.f \
  qinit.f \
  rpn3bu.f \
  rpt3bu.f \
  rptt3bu.f \
  setprob.f \
  $(CLAW)/clawpack/3d/lib/setaux.f \
  $(CLAW)/clawpack/3d/lib/claw3ez.f \
  $(CLAW)/clawpack/3d/lib/bc3.f \
  $(CLAW)/clawpack/3d/lib/b4step3.f \
  $(CLAW)/clawpack/3d/lib/chkmth.f \
  $(CLAW)/clawpack/3d/lib/claw3.f \
  $(CLAW)/clawpack/3d/lib/step3.f \
  $(CLAW)/clawpack/3d/lib/step3ds.f \
  $(CLAW)/clawpack/3d/lib/dimsp3.f \
  $(CLAW)/clawpack/3d/lib/flux3.f \
  $(CLAW)/clawpack/3d/lib/copyq3.f \
  $(CLAW)/clawpack/3d/lib/limiter.f \
  $(CLAW)/clawpack/3d/lib/philim.f \
  $(CLAW)/clawpack/3d/lib/src3.f

OUTOBJECTS = $(CLAW)/clawpack/3d/lib/out3.o \
	$(CLAW)/clawpack/3d/lib/restart3.o
OUTSOURCES = $(CLAW)/clawpack/3d/lib/out3.f \
	$(CLAW)/clawpack/3d/lib/restart3.f

HDFOBJECTS = $(CLAW)/clawpack/3d/lib/out3_hdf.o \
	$(CLAW)/clawpack/3d/lib/restart3_hdf.o
HDFSOURCES = $(CLAW)/clawpack/3d/lib/out3_hdf.f \
	$(CLAW)/clawpack/3d/lib/restart3_hdf.f

xclaw: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclaw

xclawhdf: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) $(HDFLIBS) \
	-o xclawhdf

program : $(SOURCES) $(OUTSOURCES)
	cat $(SOURCES) $(OUTSOURCES) > claw3program.f

programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > claw3program.f

clean:
	-rm -f $(OBJECTS) xclaw xclawhdf

### DO NOT remove this line - make depends on it ###
