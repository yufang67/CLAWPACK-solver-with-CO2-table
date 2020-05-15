c
c
c =========================================================
      subroutine out2(maxmx,maxmy,meqn,mbc,mx,my,
     &                 xlower,ylower,dx,dy,q,t,iframe)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 2 dimensions to a hdf file as a scientific data set.
c     # See http://hdf.ncsa.uiuc.edu/ for more info on HDF.
c     # The results in this output file can be plotted in MATLAB
c     # using the "plotclaw2" script.
c
c     # MPI+HDF Output Routines by Peter Blossey, 2003.
c     # Based on work of Donna Calhoun & Sorin Mitran.
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'
c
      parameter   (nDim = 2)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      character*14 fname
      character*13 qname2
      character*9  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id, sds_rank
      integer    sds_dims, sds_start, sds_edges, sds_stride
      dimension  id_coords(2), mstart(2), mtotal(2), idimlength(2)
      dimension  isds_id(meqn), istatus(MPI_STATUS_SIZE)
      dimension  sds_dims(nDim), sds_start(nDim)
      dimension  sds_edges(nDim), sds_stride(nDim) 

c     # HDF: x- and y- dimensions are reversed when output to HDF file.
      dimension  qbuf(21), qout(my,mx)
c
c     # HDF: Declare external HDF functions
c
      integer  sfstart, sfcreate, sfwdata, sfselect, sfendacc, sfend
      external sfstart, sfcreate, sfwdata, sfselect, sfendacc, sfend
c
c     # HDF: Set up HDF constants
c
      integer 	DFACC_READ, DFACC_WRITE, DFACC_CREATE 
      parameter(DFACC_READ = 1, DFACC_WRITE = 2, DFACC_CREATE = 4)

      integer   DFNT_FLOAT64, DFNT_INT32
      parameter(DFNT_FLOAT64 = 6, DFNT_INT32 = 24)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
      common /mpicomm/ mpi_comm_2d, lx, ly, mtotal, mstart
      common /mpi_proc_info/ np, id
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file
c
      fname = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
     &     // '.hdf'
c
c     # Specify grid number and create a string which will describe this
c     # grid in HDF file.  This could be an input for simulations with
c     # multiple grids, as in AMRCLAW.  
c
      ngrids_out=1
      qname = 'grid_'
     &     // char(ichar('0') + mod(ngrids_out/1000,10)) 
     &     // char(ichar('0') + mod(ngrids_out/100,10)) 
     &     // char(ichar('0') + mod(ngrids_out/10,10)) 
     &     // char(ichar('0') + mod(ngrids_out,10))
c
c     # MPI: Let processor which has mxstart=mystart=0 create the 
c     #      HDF file and write the grid parameters to the HDF file.
c     #      Note that this processor holds the correct global values
c     #      for xlower and ylower.
c
      if (mstart(1)+mstart(2).eq.0) then
c
c     # HDF: create hdf file.
c
         sd_id = sfstart(fname, DFACC_CREATE)
         if (sd_id.eq.FAIL) THEN
            WRITE(*,*) 
     &           'Failed to create HDF file (call to sfstart)'
            STOP
         end if
c
c     # HDF: create a data set for parameters describing q in HDF file.
c
         irank = 1
         idims = 21
      
         sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,irank,idims)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 
     &           'Failed to create scientific data set in HDF file'
            STOP
         end if
c
c     # HDF: set up parameters describing data set.
c
         istart  = 0
         iedges  = idims
         istride = 1        
         qbuf(1) = ngrids_out
         qbuf(2) = nDim
         qbuf(3) = t
         qbuf(4) = meqn
         qbuf(5) = 1.
         qbuf(6) = mtotal(1)
         qbuf(7) = mtotal(2)
         qbuf(8) = 0.
         qbuf(9) = 0.
         qbuf(10) = xlower
         qbuf(11) = ylower
         qbuf(12) = 0.
         qbuf(13) = 0.
         qbuf(14) = xlower+mtotal(1)*dx
         qbuf(15) = ylower+mtotal(2)*dy
         qbuf(16) = 0.
         qbuf(17) = 0.
         qbuf(18) = dx
         qbuf(19) = dy
         qbuf(20) = 0.
         qbuf(21) = 0.
         istat = sfwdata(sds_id,istart,istride,iedges,qbuf)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to write grid data (call to sfwdata)'
            STOP
         end if
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close grid data (call to sfendacc)'
            STOP
         end if
c
c     # HDF: Close HDF file.
c     
         istat = sfend(sd_id)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close SDS (call to sfendacc)'
            STOP
         end if
c
      end if
c
c     # MPI: Wait for grid information to be written.
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
c     # MPI: Wait for message from node 0 before writing data to HDF file.
c
      if (id.ne.0) then
         call mpi_recv(ii, 1, MPI_INTEGER, 0, id, 
     &        MPI_COMM_WORLD, istatus, ierr)
      end if
c
c     # HDF: open hdf file.
c
      sd_id = sfstart(fname, DFACC_WRITE)
      if (sd_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to open HDF file (call to sfstart)'
         STOP
      end if
c     
c     # Loop over fields in q
c     
      do m = 1,meqn
c     
c     # HDF: create a data set for parameters describing q in HDF file.
c     
         qname2 = qname // '_'
     &        // char(ichar('0') + mod(m/100,10)) 
     &        // char(ichar('0') + mod(m/10,10)) 
     &        // char(ichar('0') + mod(m,10))
c     
         if (id.eq.0) then
c
c     # HDF: Reverse dimensions when storing arrays because of different
c     #      conventions between c and FORTRAN as to which dimension should
c     #      be stored first.  Reversing the dimensions here will make
c     #      the x-direction first when reading into MATLAB.
c
            sds_rank = nDim
            sds_dims(1) = mtotal(2)
            sds_dims(2) = mtotal(1)
c     
            sds_id = sfcreate(sd_id, qname2, DFNT_FLOAT64, sds_rank,
     &           sds_dims)
            if (sds_id.eq.FAIL) THEN
               WRITE(*,*) 'Processor ', id, 
     &           '  Failed to create data set in HDF file'
               STOP
            end if
         else
c
c     # MPI/HDF: Select current data set.
c
            sds_id = sfselect(sd_id,(ngrids_out-1)*(meqn+1)+m)
            if (sds_id.eq.FAIL) THEN
               WRITE(*,*) 'Processor ', id, 
     &           '  Failed to select data set in HDF file'
               STOP
            end if
         end if
c     
c     # HDF: set up parameters describing data set.
c     
         sds_start(1)  = mstart(2)
         sds_start(2)  = mstart(1)

         sds_edges(1)  = my
         sds_edges(2)  = mx

         sds_stride(1) = 1        
         sds_stride(2) = 1        
c     
c     # Copy current field of q into qout.
c     
         do j = 1,mx
            do i = 1,my
               qout(i,j) = q(j,i,m)
            end do
         end do
c     
c     # HDF: write current processor's slab of data to hdf file.
c     
         istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qout)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Processor ', id, 
     &           '  Failed to write SDS (call to sfwdata)'
            STOP
         end if
c     
c     # HDF: Close the data set
c     
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close SDS (call to sfendacc)'
            STOP
         end if
      end do
c     
c     # HDF: Close HDF file.
c     
      istat = sfend(sd_id)
      if (istat.eq.FAIL) then
         WRITE(*,*) 'Failed to close SDS (call to sfendacc)'
         STOP
      end if
c     
c     # MPI: Master node 0 tells each processor to write to data file.
c
      if (id.eq.0) then
         do i = 1,np-1
            ii = i
c
c     # MPI: Send message to processor i, telling it to start writing.
            call mpi_send(ii, 1, MPI_INTEGER, i, i,
     &           MPI_COMM_WORLD, ierr)
c
c     # MPI: Wait for message from processor i, saying it's finished writing.
            call mpi_recv(ii, 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, 
     &           istatus, ierr)
         end do
      else
         ii = 0
c
c     # MPI: Send message to processor 0, telling it that writing is complete.
         call mpi_send(ii, 1, MPI_INTEGER, 0, id, MPI_COMM_WORLD, ierr)
      end if
      
      return
      end
