c
c
c
c     =================================================================
      subroutine claw2ez_mpi(mx,my,meqn,mbc,mwaves,maux,method,
     &      mthlim,dtv,cflv,nv,nout,outstyle,tfinal,tout,
     &      nstepout,t0,xlower,xupper,ylower,yupper,mthbc,mwork)

c     =================================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c     Documentation is available at
c                 http://www.amath.washington.edu/~claw/doc.html
c
c     Author: Peter Blossey, based on work of Donna Calhoun & Sorin Mitran
c     Version of July 2003 --  CLAWPACK Version 4.2
c
c     Author: Randall J. LeVeque
c     Version of August, 2002 --  CLAWPACK Version 4.1
c
c
      implicit double precision (a-h,o-z)

      external bc2,rpn2,rpt2,src2,b4step2

c     # MPI : These arrays are automatically allocated here.  They will be
c     # automatically deallocated when we exit this routine.
      dimension    q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      dimension  aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
      dimension work(mwork)

c     # MPI : These arrays allocated in claw2ez_mpi_driver, so they are only
c     # dimensioned here
      dimension mthlim(mwaves)
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(6)
      dimension tout(100)  !! dimensioned in claw2ez_mpi_driver.
      logical outt0
      common /restrt_block/ tinitial, iframe, outt0
      common /mpi_proc_info/ np, id

c
      if (id.eq.0) then
         open(10,file='fort.info',status='unknown',form='formatted')
      end if
c
c
c     # -----------------------------------------------------------------
c     # MPI : We have removed all code that reads in data, and does error
c     # checking on input data.  We can start right in with other stuff..
c     # -----------------------------------------------------------------

c     # MPI : Set these here so that we are sure that they are defined
      maxmx = mx
      maxmy = my

c
c     # grid spacing
      dx = (xupper - xlower) / float(mx)
      dy = (yupper - ylower) / float(my)
c


c     # time increments between outputing solution:
      if (outstyle .eq. 1) then
         dtout = (tfinal - t0)/float(nout)
      endif

c
c     # MPI : The file fort.nplot (unit=11) was created in claw2ez_mpi_driver,
C     # so we comment out the code here.

c       write(11,1101) nout
c       write(11,1101) 1
c c
c  1101 format(i5)

c     # -------------------------------------------------------------------
c     # MPI : All code from this point is taken directly from claw2ez, and
C     # has not been changed.

c
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      call setprob
c
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &        maux,aux)
c$$$         call bc2_aux(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
c$$$     &               dx,dy,q,maux,aux,told,dt,mthbc)
      endif
c
c     # set initial conditions:
c
      outt0  = .true.
      iframe = 0
      tinitial = t0
c
      call qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &     q,maux,aux)
c
c     # Reset initial time if changed within qinit (by restart2_mpi_hdf.f)
c
      t0 = tinitial
c
      outt0 = .true.
      if (outt0) then
c        # output initial data
         call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &        dx,dy,q,t0,iframe)
         if (id.eq.0) write(6,601) iframe, t0
      endif

c
c     ----------
c     Main loop:
c     ----------
c

      tend = t0
      n0   = iframe*nstepout + 1
      do 100 n=n0,nout
         tstart = tend
         if (outstyle .eq. 1)  tend = tstart + dtout
         if (outstyle .eq. 2)  tend = tout(n)
         if (outstyle .eq. 3)  tend = tstart - 1.d0  !# single-step mode
c
         call claw2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc2,rpn2,rpt2,src2,b4step2)
c
c        # check to see if an error occured:
c        if (info .ne. 0) then
c           write(6,*) '*** ERROR in claw2 ***  info =',info
c           go to 999
c           endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
c        # if outstyle=1 or 2, then nstepout=1 and we output every time
c        # we reach this point, since claw1 was called for the entire time
c        # increment between outputs.
c
c        # if outstyle=3 then we only output if we have taken nstepout
c        # time steps since the last output.

c        # iframe is the frame number used to form file names in out1
         iframe = n/nstepout
         if (iframe*nstepout .eq. n) then
            call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &            dx,dy,q,tend,iframe)
c           # only write out time and run info from master node.
            if (id.eq.0) then
               write(6,601) iframe,tend
               write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &              cflv(3),cflv(4),nv(2)
            end if
         endif
c
c        # formats for writing out information about this call to claw:

  601    format('CLAW2EZ: Frame ',i4,
     &         ' matlab plot files done at time t = ', d12.4)
c
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
  100    continue
c
  999 continue
c
      if (id.eq.0) close(10)
c
      return
      end
