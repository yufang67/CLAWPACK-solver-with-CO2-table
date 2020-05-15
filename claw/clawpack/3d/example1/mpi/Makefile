#
#               Makefile for the mpiclaw code:
#
#       To make an MPI executable, type:          make xclawmpi
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawmpihdf
#       (from the application sub-directory)
#
#       To combine ASCII output files (one from
#       each processor at each time level) into
#       fort.qXXXX files for use with MATLAB, type:  make catfiles
#       (in directory with fort.qXXXX.YY files)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
#       To combine all source files including
#       HDF output routines type:                 make programhdf
#       (from the application sub-directory)
#
FFLAGS = -g -O3
LFLAGS = $(FFLAGS)
F77    = mpif77 -c -mismatch # -mismatch is necessary for NAG f95 compiler.
LINK   = mpif77
HDFLIBS = -L/usr/local/hdf/lib -lmfhdf -ldf -ljpeg -lz

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  ../qinit.o \
  ../rpn3bu.o \
  ../rpt3bu.o \
  ../rptt3bu.o \
  ../setprob.o

LIBOBJECTS = \
  $(CLAW)/clawpack/3d/mpilib/claw3ez_mpi_driver.o \
  $(CLAW)/clawpack/3d/mpilib/claw3ez_mpi.o \
  $(CLAW)/clawpack/3d/mpilib/bc3_mpi.o \
  $(CLAW)/clawpack/3d/mpilib/bc3_aux_mpi.o \
  $(CLAW)/clawpack/3d/mpilib/claw3_mpi.o \
  $(CLAW)/clawpack/3d/lib/setaux.o \
  $(CLAW)/clawpack/3d/lib/b4step3.o \
  $(CLAW)/clawpack/3d/lib/chkmth.o \
  $(CLAW)/clawpack/3d/lib/step3.o \
  $(CLAW)/clawpack/3d/lib/step3ds.o \
  $(CLAW)/clawpack/3d/lib/dimsp3.o \
  $(CLAW)/clawpack/3d/lib/flux3.o \
  $(CLAW)/clawpack/3d/lib/copyq3.o \
  $(CLAW)/clawpack/3d/lib/limiter.o \
  $(CLAW)/clawpack/3d/lib/philim.o \
  $(CLAW)/clawpack/3d/lib/src3.o

SOURCES = \
  ../qinit.f \
  ../rpn3bu.f \
  ../rpt3bu.f \
  ../rptt3bu.f \
  ../setprob.f \
  $(CLAW)/clawpack/3d/mpilib/claw3ez_mpi_driver.f \
  $(CLAW)/clawpack/3d/mpilib/claw3ez_mpi.f \
  $(CLAW)/clawpack/3d/mpilib/bc3_mpi.f \
  $(CLAW)/clawpack/3d/mpilib/bc3_aux_mpi.f \
  $(CLAW)/clawpack/3d/mpilib/claw3_mpi.f \
  $(CLAW)/clawpack/3d/lib/setaux.f \
  $(CLAW)/clawpack/3d/lib/b4step3.f \
  $(CLAW)/clawpack/3d/lib/chkmth.f \
  $(CLAW)/clawpack/3d/lib/step3.f \
  $(CLAW)/clawpack/3d/lib/step3ds.f \
  $(CLAW)/clawpack/3d/lib/dimsp3.f \
  $(CLAW)/clawpack/3d/lib/flux3.f \
  $(CLAW)/clawpack/3d/lib/copyq3.f \
  $(CLAW)/clawpack/3d/lib/limiter.f \
  $(CLAW)/clawpack/3d/lib/philim.f \
  $(CLAW)/clawpack/3d/lib/src3.f

OUTOBJECTS = $(CLAW)/clawpack/3d/mpilib/out3_mpi.o \
	$(CLAW)/clawpack/3d/mpilib/restart3_mpi.o
OUTSOURCES = $(CLAW)/clawpack/3d/mpilib/out3_mpi.f \
	$(CLAW)/clawpack/3d/mpilib/restart3_mpi.f

HDFOBJECTS = $(CLAW)/clawpack/3d/mpilib/out3_mpi_hdf.o \
	$(CLAW)/clawpack/3d/mpilib/restart3_mpi_hdf.o
HDFSOURCES = $(CLAW)/clawpack/3d/mpilib/out3_mpi_hdf.f \
	$(CLAW)/clawpack/3d/mpilib/restart3_mpi_hdf.f

xclawmpi: $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclawmpi 

xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xclawmpi \
	$(HDFLIBS)

program : $(SOURCES) $(OUTSOURCES)
	cat $(SOURCES) $(OUTSOURCES) > claw3program.f

programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > claw3program.f

# Type "make catfiles" to combine fort.qXXXX.YY output files from out2_mpi
# into fort.qXXXX files (for use with CLAWPACK's MATLAB graphics routines).
# This step is not necessary with the HDF output routines.
catfiles : fort.q0000

fort.q0000 : fort.q0000.00
	$(CLAW)/clawpack/3d/mpilib/catfiles

fort.q0000.00 :
	$(CLAW)/clawpack/3d/mpilib/catfiles

### DO NOT remove this line - make depends on it ###
