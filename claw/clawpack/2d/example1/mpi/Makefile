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
FFLAGS = -O3
LFLAGS = $(FFLAGS)
F77    = mpif77 -c # -mismatch is necessary for NAG f95 compiler.
LINK   = mpif77
HDFLIBS = -L/usr/local/hdf/lib -lmfhdf -ldf -ljpeg -lz

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  ../rpn2.o \
  ../rpt2.o \
  ../qinit.o \
  ../setprob.o

LIBOBJECTS = \
  $(CLAW)/clawpack/2d/mpilib/claw2ez_mpi_driver.o \
  $(CLAW)/clawpack/2d/mpilib/claw2ez_mpi.o \
  $(CLAW)/clawpack/2d/mpilib/claw2_mpi.o \
  $(CLAW)/clawpack/2d/mpilib/bc2_mpi.o \
  $(CLAW)/clawpack/2d/mpilib/bc2_aux_mpi.o \
  $(CLAW)/clawpack/2d/lib/setaux.o \
  $(CLAW)/clawpack/2d/lib/b4step2.o \
  $(CLAW)/clawpack/2d/lib/step2.o \
  $(CLAW)/clawpack/2d/lib/step2ds.o \
  $(CLAW)/clawpack/2d/lib/dimsp2.o \
  $(CLAW)/clawpack/2d/lib/flux2.o \
  $(CLAW)/clawpack/2d/lib/copyq2.o \
  $(CLAW)/clawpack/2d/lib/limiter.o \
  $(CLAW)/clawpack/2d/lib/philim.o \
  $(CLAW)/clawpack/2d/lib/src2.o  
  
SOURCES = \
  ../rpn2.f \
  ../rpt2.f \
  ../qinit.f \
  ../setprob.f \
  $(CLAW)/clawpack/2d/mpilib/claw2ez_mpi_driver.f \
  $(CLAW)/clawpack/2d/mpilib/claw2ez_mpi.f \
  $(CLAW)/clawpack/2d/mpilib/claw2_mpi.f \
  $(CLAW)/clawpack/2d/mpilib/bc2_mpi.f \
  $(CLAW)/clawpack/2d/mpilib/bc2_aux_mpi.f \
  $(CLAW)/clawpack/2d/lib/setaux.f \
  $(CLAW)/clawpack/2d/lib/b4step2.f \
  $(CLAW)/clawpack/2d/lib/step2.f \
  $(CLAW)/clawpack/2d/lib/step2ds.f \
  $(CLAW)/clawpack/2d/lib/dimsp2.f \
  $(CLAW)/clawpack/2d/lib/flux2.f \
  $(CLAW)/clawpack/2d/lib/copyq2.f \
  $(CLAW)/clawpack/2d/lib/limiter.f \
  $(CLAW)/clawpack/2d/lib/philim.f \
  $(CLAW)/clawpack/2d/lib/src2.f  
  
OUTOBJECTS = $(CLAW)/clawpack/2d/mpilib/out2_mpi.o \
	$(CLAW)/clawpack/2d/mpilib/restart2_mpi.o
OUTSOURCES = $(CLAW)/clawpack/2d/mpilib/out2_mpi.f \
	$(CLAW)/clawpack/2d/mpilib/restart2_mpi.f

HDFOBJECTS = $(CLAW)/clawpack/2d/mpilib/out2_mpi_hdf.o \
	$(CLAW)/clawpack/2d/mpilib/restart2_mpi_hdf.o
HDFSOURCES = $(CLAW)/clawpack/2d/mpilib/out2_mpi_hdf.f \
	$(CLAW)/clawpack/2d/mpilib/restart2_mpi_hdf.f

xclawmpi: $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclawmpi 

xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xclawmpi \
	$(HDFLIBS)

program : $(SOURCES) $(OUTSOURCES)
	cat $(SOURCES) $(OUTSOURCES) > claw2program.f

programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > claw2program.f

# Type "make catfiles" to combine fort.qXXXX.YY output files from out2_mpi
# into fort.qXXXX files (for use with CLAWPACK's MATLAB graphics routines).
# This step is not necessary with the HDF output routines.
catfiles : fort.q0000

fort.q0000 : fort.q0000.00
	$(CLAW)/clawpack/3d/mpilib/catfiles

fort.q0000.00 :
	$(CLAW)/clawpack/3d/mpilib/catfiles

### DO NOT remove this line - make depends on it ###
