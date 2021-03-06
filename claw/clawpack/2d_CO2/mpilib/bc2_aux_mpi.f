
c
c
c     =====================================================
      subroutine bc2_aux(maxmx,maxmy,meqn,mbc,mx,my,xlower,
     &               ylower,dx,dy,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (xlower),  2 (xupper),
c     #                       3 (ylower),  4 (yupper):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to a layer of mbc ghost cells outside the region.
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'

      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

c     # MPI : Create arrays for communication with neighboring processors
c     # These will be automatically allocated here, and will be deallocated
c     # when we exit this routine.

      dimension ql(mbc, 1:my, maux)
      dimension qr(mbc, 1:my, maux)
      dimension qt(1-mbc:mx+mbc, mbc, maux)
      dimension qb(1-mbc:mx+mbc, mbc, maux)

      dimension ql_send(mbc, 1:my, maux)
      dimension qr_send(mbc, 1:my, maux)
      dimension qt_send(1-mbc:mx+mbc, mbc, maux)
      dimension qb_send(1-mbc:mx+mbc, mbc, maux)

      parameter (ileft = 1, iright = 2, ibottom = 3, itop = 4)

      dimension mthbc(4), ireq(4), MPIstatus(MPI_STATUS_SIZE,4)
      dimension id_coords(2), id_next(2), id_last(2)
      dimension mtotal(2), mstart(2)

c     # mpi : this communicator gives us information about the mapping 
c     # of nodes to grid subdomains.
      common /mpicomm/ mpi_comm_2d, lx, ly, mtotal, mstart
      common /mpi_proc_info/ np, id
c
c
c     # mpi : initialize number of outstanding communication requests
      nreq = 0

c     # Get the coordinates of the current processor in the processor array.
c     # id_coords(1) is the i coordinate of this processor in the array.
c     # id_coords(2) is the j coordinate of this processor in the array.
      call mpi_cart_coords(mpi_comm_2d,id,2,id_coords,ierr)
      imap = id_coords(1)
      jmap = id_coords(2)
c
c-------------------------------------------------------
c     # left boundary (xlower):
c-------------------------------------------------------
c     # MPI : added option 4 : internal boundary
      go to (100,110,120,130,140) mthbc(ileft) + 1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      go to 199
c
  110 continue
c     # zero-order extrapolation: nothing to do for aux arrays
      go to 199

  120 continue
c     # MPI : periodic
      if (lx == 1) then
         do m = 1,maux
            do j = 1,my
               do i = 1,mbc
                  ql(i,j,m) = aux(mx-mbc+i,j,m)
               end do
            end do
         end do
      else
         STOP 'periodic boundaries are treated as internal for lx>1'
      endif
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
c     # negate the normal velocity:
      go to 199

  140 continue
c     # MPI : Internal boundary

c     # Send left side of domain
      nbcdata = mbc*my*maux
      id_last(1) = imap-1
      id_last(2) = jmap
      call mpi_cart_rank(mpi_comm_2d,id_last,idcomm,ierr)

c     # Receive left ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq + 1
      call mpi_irecv(ql, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, ileft, MPI_COMM_WORLD, ireq(nreq), ierr)

      do m = 1,maux
         do j = 1,my
            do i = 1,mbc
               ql_send(i,j,m) = aux(i,j,m)
            end do
         end do
      end do
      call mpi_send(ql_send, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, iright, MPI_COMM_WORLD, ierr)

      go to 199
  199 continue
c
c-------------------------------------------------------
c     # right boundary (xupper):
c-------------------------------------------------------
c     # MPI : added option 4 : internal boundary
      go to (200,210,220,230,240) mthbc(iright)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      go to 299

  210 continue
c     # zero-order extrapolation:
      go to 299

  220 continue
c     # MPI : periodic:
      if (lx == 1) THEN
c        # All data is local. Apply periodic boundary condition.
         do m = 1,maux
            do j = 1,my
               do i = 1,mbc
                  qr(i,j,m) = aux(i,j,m)
               end do
            end do
         end do
      else
         STOP 'periodic boundaries are treated as internal for lx>1'
      endif
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
c     # negate the normal velocity:
      go to 299

  240 continue
c     # MPI : Internal boundary condition

c     # Send right side of domain
      nbcdata = mbc*my*maux
      id_next(1) = imap+1
      id_next(2) = jmap
      call mpi_cart_rank(mpi_comm_2d,id_next,idcomm,ierr)

c     # Receive right ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq+1
      call mpi_irecv(qr, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, iright, MPI_COMM_WORLD, ireq(nreq), ierr)

      do m = 1,maux
         do j = 1,my
            do i = 1,mbc
               qr_send(i,j,m) = aux(mx-mbc+i,j,m)
            end do
         end do
      end do
      call mpi_send(qr_send,nbcdata,MPI_DOUBLE_PRECISION,
     &      idcomm, ileft, MPI_COMM_WORLD, ierr)

      go to 299
  299 continue
c
c-------------------------------------------------------
c     # MPI: If information has been passed from other processors, 
c     # record it in q so that corner ghost cells in y- and z-
c     # directions will have correct information.
c-------------------------------------------------------
c     # MPI : Wait for all processes to complete
      do i = 1,nreq
         call mpi_wait(ireq(i),MPIstatus(1,i),ierr)
      enddo
c     # MPI : re-initialize number of outstanding communication requests
      nreq = 0

c     # Copy values from communication buffers back into main array
c     # Left boundary
      select case (mthbc(ileft))
      case (2,4)
         do m = 1,maux
            do j = 1,my
               do i = 1,mbc
                  aux(i-mbc,j,m) = ql(i,j,m)
               end do
            end do
         end do
      case default
      end select

c     # Right boundary
      select case (mthbc(iright))
      case (2,4)
         do m = 1,maux
            do j = 1,my
               do i = 1,mbc
                  aux(mx+i,j,m) = qr(i,j,m)
               end do
            end do
         end do
      case default
      end select
c
c-------------------------------------------------------
c     # bottom boundary (ylower):
c-------------------------------------------------------
c     # MPI : added option 4 : internal boundary
      go to (300,310,320,330,340) mthbc(ibottom)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      go to 399

  320 continue
c     # MPI : periodic:
      if (ly == 1) THEN
c        # All data is local. Apply periodic boundary condition.
         do m = 1,maux
            do j = 1,mbc
               do i = 1-mbc,mx+mbc
                  qb(i,j,m) = aux(i,my-mbc+j,m)
               end do
            end do
         end do
      else
         STOP 'periodic boundaries are treated as internal for ly>1'
      endif
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
c     # negate the normal velocity:
      go to 399

  340 continue
c     # MPI : Internal boundary condition

c     # Send bottom side of domain
      nbcdata = (mx+2*mbc)*mbc*maux
      id_last(1) = imap
      id_last(2) = jmap-1
      call mpi_cart_rank(mpi_comm_2d,id_last,idcomm,ierr)

c     # Receive bottom ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq+1
      CALL mpi_irecv(qb, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, ibottom, MPI_COMM_WORLD, ireq(nreq), ierr)

      do m = 1,maux
         do j = 1,mbc
            do i = 1-mbc,mx+mbc
               qb_send(i,j,m) = aux(i,j,m)
            end do
         end do
      end do
      call mpi_send(qb_send, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm,itop, MPI_COMM_WORLD, ierr)

      go to 399
  399 continue
c
c-------------------------------------------------------
c     # top boundary (yupper):
c-------------------------------------------------------
c     # MPI : added option 4 : internal boundary
      go to (400,410,420,430,440) mthbc(itop)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      go to 499

  410 continue
c     # zero-order extrapolation:
      go to 499

  420 continue
c     # MPI : periodic:
      if(ly == 1) THEN
c        # All data is local. Apply periodic boundary condition.
         do m = 1,maux
            do j = 1,mbc
               do i = 1-mbc,mx+mbc
                  qt(i,j,m) = aux(i,j,m)
               end do
            end do
         end do
      else
         STOP 'periodic boundaries are treated as internal for ly>1'
      endif
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
c     # negate the normal velocity:
      go to 499

  440 continue
c     # MPI : Internal boundary condition

c     # Send top of domain
      nbcdata = (mx+2*mbc)*mbc*maux
      id_next(1) = imap
      id_next(2) = jmap+1
      call mpi_cart_rank(mpi_comm_2d,id_next,idcomm,ierr)

c     # Receive top ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq+1
      call mpi_irecv(qt, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm,itop, MPI_COMM_WORLD, ireq(nreq), ierr)

      do m = 1,maux
         do j = 1,mbc
            do i = 1-mbc,mx+mbc
               qt_send(i,j,m) = aux(i,my-mbc+j,m)
            end do
         end do
      end do
      call mpi_send(qt_send,nbcdata,MPI_DOUBLE_PRECISION,
     &      idcomm,ibottom, MPI_COMM_WORLD, ierr)

      go to 499
  499 continue

c
c-------------------------------------------------------
c     # MPI: If information has been passed from other processors, 
c     # record it in q so that corner ghost cells in z-direction
c     # will have correct information.
c-------------------------------------------------------
c     # MPI : Wait for all processes to complete
      do i = 1,nreq
         call mpi_wait(ireq(i),MPIstatus(1,i),ierr)
      enddo
c     # MPI : re-initialize number of outstanding communication requests
      nreq = 0

c     # Copy values from communication buffers back into main array
c     #   Bottom boundary
      select case (mthbc(ibottom))
      case (2,4)
         do m = 1,maux
            do j = 1,mbc
               do i = 1-mbc,mx+mbc
                  aux(i,j-mbc,m) = qb(i,j,m)
               end do
            end do
         end do
      case default
      end select

c     # Top boundary
      select case (mthbc(itop))
      case (2,4)
         do m = 1,maux
            do j = 1,mbc
               do i = 1-mbc,mx+mbc
                  aux(i,my+j,m) = qt(i,j,m)
               end do
            end do
         end do
      case default
      end select

      return
      end
