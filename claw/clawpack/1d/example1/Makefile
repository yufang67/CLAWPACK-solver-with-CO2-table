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
#       To compile the lib library files type:    make 
#       (from the claw/clawpack/1d/lib sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
FFLAGS = -O
LFLAGS = $(FFLAGS)
F77    = gfortran -c
LINK   = gfortran
HDFLIBS = -L/usr/local/hdf/lib -lmfhdf -ldf -ljpeg -lz

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  driver.o \
  qinit.o \
  rp1.o \
  setprob.o 

LIBOBJECTS = \
  $(CLAW)/clawpack/1d/lib/claw1ez.o \
  $(CLAW)/clawpack/1d/lib/bc1.o \
  $(CLAW)/clawpack/1d/lib/setaux.o \
  $(CLAW)/clawpack/1d/lib/b4step1.o \
  $(CLAW)/clawpack/1d/lib/claw1.o \
  $(CLAW)/clawpack/1d/lib/step1.o \
  $(CLAW)/clawpack/1d/lib/copyq1.o \
  $(CLAW)/clawpack/1d/lib/limiter.o \
  $(CLAW)/clawpack/1d/lib/philim.o \
  $(CLAW)/clawpack/1d/lib/src1.o 

SOURCES = \
  driver.f \
  qinit.f \
  rp1.f \
  setprob.f \
  $(CLAW)/clawpack/1d/lib/claw1ez.f \
  $(CLAW)/clawpack/1d/lib/bc1.f \
  $(CLAW)/clawpack/1d/lib/setaux.f \
  $(CLAW)/clawpack/1d/lib/b4step1.f \
  $(CLAW)/clawpack/1d/lib/claw1.f \
  $(CLAW)/clawpack/1d/lib/step1.f \
  $(CLAW)/clawpack/1d/lib/copyq1.f \
  $(CLAW)/clawpack/1d/lib/limiter.f \
  $(CLAW)/clawpack/1d/lib/philim.f \
  $(CLAW)/clawpack/1d/lib/src1.f 

OUTOBJECTS = $(CLAW)/clawpack/1d/lib/out1.o
OUTSOURCES = $(CLAW)/clawpack/1d/lib/out1.f

HDFOBJECTS = $(CLAW)/clawpack/1d/lib/out1_hdf.o
HDFSOURCES = $(CLAW)/clawpack/1d/lib/out1_hdf.f

xclaw: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclaw

xclawhdf: $(OBJECTS) 
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) $(HDFLIBS) \
	-o xclawhdf

program : $(SOURCES) $(OUTSOURCES)
	cat $(SOURCES) $(OUTSOURCES) > claw1program.f

programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > claw1program.f

clean:
	-rm -f $(OBJECTS) xclawhdf

### DO NOT remove this line - make depends on it ###
