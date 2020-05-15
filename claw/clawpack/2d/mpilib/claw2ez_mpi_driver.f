c     
c     
c     
c     =================================================================
      program claw2ez_mpi_driver
c     =================================================================
c     
c     Author: Peter Blossey, based on work of Donna Calhoun & Sorin Mitran
c     Version of July 2003 --  CLAWPACK Version 4.2
c
c     This is the "driver" routine for claw2ez_mpi, which makes the calls to
c     claw2, etc.  This file reads in claw2ez.data, gets all the relevant
c     dimensions, and then passes them off to claw2ez_mpi.f, where memory for
c     q, aux and work is automatically allocated.
c     
c     
      implicit double precision (a-h,o-z)
      include 'mpif.h'

      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(6),mthbcp(6)
      dimension tout(100), mthlim(100)
      dimension int_param(30), dble_param(30)
      dimension ndims(2), id_coord(2), mtotal(2), mstart(2)
      logical   periods(2), reorder

c     # MPI : This common block is used by out2_mpi and bc2_mpi
      common /mpicomm/ mpi_comm_2d, lx, ly, mtotal, mstart
      common /mpi_proc_info/ np, id
c
c     # ------------------------------------
c     # MPI : Calls to set up MPI process go here.
      call mpi_init(ierr)
      if (ierr .ne. mpi_success) then
         write(6,*) 'Unable to initialize MPI system'
         call mpi_finalize(ierr)
         stop
      else
         call mpi_comm_size(mpi_comm_world,np,ierr)
         call mpi_comm_rank(mpi_comm_world,id,ierr)
      endif

      if (id.eq.0) then

c     # MPI : Read in data on master node before passing to other processors
         open(55,file='claw2ez.data',status='old',form='formatted')

c     # ----------------------------------------------------------
c     # MPI : Read in data from claw2ez.data.  This part of the code is
c     #       essentially a cut and paste from claw2ez.f, with the following
c     #       differences :
c     #
c     #        - meqn1, mwaves1 and mbc1 are changed to meqn, mwaves and mbc
c     #        - we set maux = method(7)
c     #        - checks involving meqn1, mwaves1, mbc1, maux have been removed
c     #        - all references to maxmx, maxmy and maxmz have been changed to
C     #          mx, my and mz (or mxp, myp, mzp)
c     #        - Code has been added to divide up domain among processors
c     # ----------------------------------------------------------


c     # Read the input in standard form from claw2ez.data:

c     domain variables
         read(55,*) mx
         read(55,*) my

c     i/o variables
         read(55,*) nout
         read(55,*) outstyle
         if (outstyle.eq.1) then
            read(55,*) tfinal
            nstepout = 1
         elseif (outstyle.eq.2) then
            read(55,*) (tout(i), i=1,nout)
            nstepout = 1
         elseif (outstyle.eq.3) then
            read(55,*) nstepout, nstop
            nout = nstop
         endif


c     timestepping variables
         read(55,*) dtv(1)
         read(55,*) dtv(2)
         read(55,*) cflv(1)
         read(55,*) cflv(2)
         read(55,*) nv(1)
c     

c     # input parameters for clawpack routines
         read(55,*) method(1)
         read(55,*) method(2)
         read(55,*) method(3)
         read(55,*) method(4)
         read(55,*) method(5)
         read(55,*) method(6)
         read(55,*) method(7)
c     # MPI : set maux directly.
         maux = method(7)

c     # MPI : set values of meqn, mwaves directly, rather than reading in
C     # meqn1, mwaves1
         read(55,*) meqn
         read(55,*) mwaves
         read(55,*) (mthlim(mw), mw=1,mwaves)

         read(55,*) t0
         read(55,*) xlower
         read(55,*) xupper
         read(55,*) ylower
         read(55,*) yupper
c     
c     # MPI : Read in mbc directly, rather than mbc1
         read(55,*) mbc
         read(55,*) mthbc(1)
         read(55,*) mthbc(2)
         read(55,*) mthbc(3)
         read(55,*) mthbc(4)

c     # MPI : Read in number of processors in each direction
         read(55,*) lx
         read(55,*) ly

c     # MPI : Done reading in data.
         close(55)

c     # ----------------------------------------
c     # MPI : Check input data

c     # Check for consistency of method specification.
c     # Note that info is initialized to zero and is reset
c     # if the method specification is inconsistent.
         info = 0

         if (mbc.lt.2) then
            write(6,*) '*** ERROR *** mbc < 2'
            write(6,*) 'need at least two ghost cells along each edge'
            info = 6
         end if  
         
         if (method(1).eq.0) then
c     # Fixed size time steps.  Compute the number of steps:
            dt = dtv(1)         
            maxn = (tfinal - t0 + 1d-10) / dt
            if ((outstyle.eq.1) .and. 
     &           (dabs(maxn*dt - (tfinal-t0)).gt.1d-8)) then
c     # dt doesn't divide time interval integer number of times      
               write(6,*) '*** ERROR *** '
               write(6,*) 'dt does not divide time interval'
               info = 6
            endif
         endif
  
         if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
c     # target cfl is larger than max cfl.
            write(6,*) '*** ERROR ***  cflv(2) > cflv(1)'
            write(6,*) 'target cfl is larger than max cfl'
            info = 6
         endif
  
         if (method(6).gt.method(7)) then  
            write(6,*) '*** ERROR ***  method(6) > method(7)'
            write(6,*) '*** ERROR ***  method(6) > method(7)'
            info = 6
         endif

         if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &        (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
            write(6,*) '*** ERROR ***  periodic boundary conditions'
            write(6,*) 'require mthbc(1) and mthbc(2) BOTH be set to 2'
            info = 6
         endif

         if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &        (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
            write(6,*) '*** ERROR ***  periodic boundary conditions'
            write(6,*) 'require mthbc(3) and mthbc(4) BOTH be set to 2'
            info = 6
         endif

c     # ----------------------------------------------------------------
c     # Check that number of processors equals lx*ly
c     # ---------------------------------------------------------------
         if (lx*ly .ne. np) then
            write(6,'(A,i2,A,i2,A,i2,A,i2)')
     &           'Number of processes (np = ',np,
     &           ') should equal number of subdomains (', lx,
     &           'x',ly,' = ', lx*ly,')'
            info = 6
         endif

c     # Checks for consistency of values meqn1, mwaves1, mbc1 and maux
c     # have been deleted from this version
      end if

c     # MPI: Broadcast info variable from master node.
c     # MPI: Terminate all processes if any data check failed.
      CALL MPI_BCAST(info,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( info .eq. 6) then
         call mpi_finalize(ierr)
         stop
      endif

      if (id.eq. 0) then
c     # MPI:  Fill two separate arrays with double precission
c     #       and integer parameters.
         int_param(1) = mx
         int_param(2) = my
         int_param(3) = 0
         int_param(4) = nout
         int_param(5) = nv(1)
         int_param(6) = method(1)
         int_param(7) = method(2)
         int_param(8) = method(3)
         int_param(9) = method(4)
         int_param(10)= method(5)
         int_param(11)= method(6)
         int_param(12)= method(7)
         int_param(13)= meqn
         int_param(14)= mwaves
         int_param(15)= mbc
         int_param(16)= mthbc(1)
         int_param(17)= mthbc(2)
         int_param(18)= mthbc(3)
         int_param(19)= mthbc(4)
         int_param(20)= 0
         int_param(21)= 0
         int_param(22)= lx
         int_param(23)= ly
         int_param(24)= 0

         
         dble_param(1) = outstyle
         dble_param(2) = dtv(1)
         dble_param(3) = dtv(2)
         dble_param(4) = cflv(1)
         dble_param(5) = cflv(2)
         dble_param(6) = t0
         dble_param(7) = xlower
         dble_param(8) = xupper
         dble_param(9) = ylower
         dble_param(10)= yupper
         dble_param(11)= 0.
         dble_param(12)= 0.
         dble_param(13)= dt
      endif

c     # MPI: Broadcast values of parameters from claw2ez.data      
      CALL MPI_BCAST(int_param,24,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dble_param,13,MPI_DOUBLE_PRECISION,
     &     0,MPI_COMM_WORLD,ierr)

      if (id.ne.0) then
c     # MPI: Set input parameters to those read in by master process
         mx = int_param(1)
         my = int_param(2)

         nout = int_param(4)
         nv(1) = int_param(5)

         method(1) = int_param(6)
         method(2) = int_param(7)
         method(3) = int_param(8)
         method(4) = int_param(9)
         method(5) = int_param(10)
         method(6) = int_param(11)
         method(7) = int_param(12)
c     # MPI : set maux directly.
         maux = method(7)

         meqn = int_param(13)
         mwaves = int_param(14)

         mbc = int_param(15)
         mthbc(1) = int_param(16)
         mthbc(2) = int_param(17)
         mthbc(3) = int_param(18)
         mthbc(4) = int_param(19)

         lx = int_param(22)
         ly = int_param(23)

         outstyle = dble_param(1)
         dtv(1)   = dble_param(2) 
         dtv(2)   = dble_param(3) 
         cflv(1)  = dble_param(4) 
         cflv(2)  = dble_param(5) 
         t0       = dble_param(6) 
         xlower   = dble_param(7) 
         xupper   = dble_param(8) 
         ylower   = dble_param(9) 
         yupper   = dble_param(10)
         dt       = dble_param(13)
      endif

c     # MPI: Broadcast output parameters as necessary.
      if (outstyle.eq.1) then
         CALL MPI_BCAST(tfinal,1,MPI_DOUBLE_PRECISION, 
     &        0,MPI_COMM_WORLD,ierr)
         nstepout = 1
      elseif (outstyle.eq.2) then
         CALL MPI_BCAST(tout,nout,MPI_DOUBLE_PRECISION,
     &        0,MPI_COMM_WORLD,ierr)
         nstepout = 1
      elseif (outstyle.eq.3) then
         CALL MPI_BCAST(nstepout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(nstop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         nout = nstop
      endif

c     # MPI: Broadcast limiter settings
      CALL MPI_BCAST(mthlim,mwaves,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
c     # ----------------------------------------------------------------
c     # Set up processor array
c     # ----------------------------------------------------------------

c     # Use MPI to set up processor map.  We want a cartesian array of
c     # processors which is lx x ly x lz.  Set up some parameters first.
c     # Specify number of processors in each dimension of array.
      ndims(1) = lx
      ndims(2) = ly
c     # Specify in which directions (if any) the array is periodic.
      do i = 1,2
         periods(i) = .false.
         if (mthbc(2*i).eq.2) periods(i) = .true.
      end do
c     # Specify whether the array can be reshuffled to locate neighboring
c     # processors on the array on the same physical machine.
      reorder = .true.
c     # Set up a new MPI commumicator -- mpi_comm_2d -- which holds info 
c     # about the array.
      call mpi_cart_create(MPI_COMM_WORLD,2,ndims,periods,reorder,
     &     mpi_comm_2d,ierr)

c     # Get the coordinates of the current processor in the processor array.
c     # id_coord(1) is the i coordinate of this processor in the array.
c     # id_coord(2) is the j coordinate of this processor in the array.
      call mpi_cart_coords(mpi_comm_2d,id,2,id_coord,ierr)

c     # ----------------------------------------------------------------
c     # Specify size of data on current processor and its global position.
c     # ----------------------------------------------------------------

      mxp = int(mx/lx)
      leftover = mx - mxp*lx
      if (id_coord(1).lt.leftover) mxp = mxp + 1
      mxstart = id_coord(1)*int(mx/lx)
      do i = 0,id_coord(1)-1
         if (i .lt. leftover) mxstart = mxstart + 1
      end do

      myp = int(my/ly)
      leftover = my - myp*ly
      if (id_coord(2).lt.leftover) myp = myp + 1
      mystart = id_coord(2)*int(my/ly)
      do i = 0,id_coord(2)-1
         if (i .lt. leftover) mystart = mystart + 1
      end do
c
c     # MPI: Save size of grid that includes all processors
c     #      as well as starting coordinates of each processor's
c     #      block of data.
c
      mtotal(1) = mx
      mtotal(2) = my

      mstart(1) = mxstart
      mstart(2) = mystart

c     # ----------------------------------------------------------------
c     # Set internal boundary conditions.
c     # ----------------------------------------------------------------

c     # MPI : Initially assume each boundary is internal.
      do i = 1,4
         mthbcp(i) = 4
      end do

c     # Use mthbc from above if we are at a physical boundary.

c     # If left or right boundary is physical, reset mthbcp
      if (id_coord(1) .eq. 0)    mthbcp(1) = mthbc(1)
      if (id_coord(1) .eq. lx-1) mthbcp(2) = mthbc(2)

c     # MPI: Treat periodic boundaries as internal.
c     #      Note that mpi_cart_rank routines supports periodic
c     #      arrays of processors if periodicity is specified in
c     #      call to mpi_cart_create.
      if ((lx.gt.1).and.(mthbc(1).eq.2)) then
         mthbcp(1) = 4
         mthbcp(2) = 4
      end if

c     # If front or back boundary is physical, reset mthbcp
      if (id_coord(2) .eq. 0)    mthbcp(3) = mthbc(3)
      if (id_coord(2) .eq. ly-1) mthbcp(4) = mthbc(4)

c     # MPI: Treat periodic boundaries as internal.
      if ((ly.gt.1).and.(mthbc(3).eq.2)) then
         mthbcp(3) = 4
         mthbcp(4) = 4
      end if

c     # -----------------------------------------------------------
c     # MPI : We set value of mwork
c     # The space for the work array will be allocated automatically in
c     # claw2ez_mpi, and so is not allocated here.

      if (method(5).lt.2) then
c        # only need one qwork array
         narray = 1
      else
c        # need two qwork arrays for Strang splitting
         narray = 2
      endif

      maxm = max0(mxp, myp)
      mwork = (maxm+2*mbc)*(10*meqn + mwaves + meqn*mwaves 
     &                      + 3*maux + 2) 
     &          + narray * (mxp + 2*mbc) * (myp + 2*mbc) * meqn   
c
c     # -------------------------------------------------------
c     # Set up upper and lower bounds for each subdomain
c     # -------------------------------------------------------
c
      xlowerp = xlower + (xupper - xlower)*dble(mxstart)/dble(mx)
      xupperp = xlower + (xupper - xlower)*dble(mxstart+mxp)/dble(mx)
      ylowerp = ylower + (yupper - ylower)*dble(mystart)/dble(my)
      yupperp = ylower + (yupper - ylower)*dble(mystart+myp)/dble(my)

c     # ---------------------------------------------------------
c     # Create file fort.nplot
c     # ---------------------------------------------------------
c     # MPI : Only node 0 needs to do this
      if  (id .eq. 0) then
         open(11,file='fort.nplot',status='unknown',form='formatted')
         write(11,'(i5)') nout
         write(11,'(i5)') 1
         close(11)
      endif

c     # ------------------------------------------------------------
c     # Call claw2ez_mpi
c     # -----------------------------------------------------------
c     # MPI : Now we send requests to each processor.  Node 0 will keep track
c     # of the time.  At this point, each processor should be able to run
c     # independently of the others.
      if (id .eq. 0) then
         cput0 = mpi_wtime()

         write(6,'(A,I5)') 'running on node ', id
         call claw2ez_mpi(mxp,myp,meqn,mbc,mwaves,maux,method,
     &         mthlim,dtv,cflv,nv,nout,outstyle,tfinal,tout,nstepout,
     &         t0,xlowerp,xupperp,ylowerp,yupperp,mthbcp,mwork)

         cput1 = mpi_wtime()
         write(6,'(A,F10.3)') 'Elapsed time = ', cput1 - cput0
      else
         write(6,'(A,I5)') 'running on node ', id
         call claw2ez_mpi(mxp,myp,meqn,mbc,mwaves,maux,method,
     &         mthlim,dtv,cflv,nv,nout,outstyle,tfinal,tout,nstepout,
     &         t0,xlowerp,xupperp,ylowerp,yupperp,mthbcp,mwork)
      endif

      call mpi_finalize(ierr)

      stop
      end
