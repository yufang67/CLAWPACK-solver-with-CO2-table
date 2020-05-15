c
c
c ===================================================================
      subroutine step1(maxmx,meqn,mwaves,mbc,mx,q,aux,dx,dt,
     &              method,mthlim,cfl,f,fwave,s,amdq,apdq,dtdx,rp1)
c ===================================================================
c
c     # Take one time step, updating q.
c
c     ----------------------------------------------------------------------
c     # step1fw is a modified version of step1 to use fwave instead of wave.
c     # A modified Riemann solver rp1 must be used in conjunction with this
c     # routine, which returns fwave's instead of wave's.
c     # See http://amath.washington.edu/~claw/fwave.html
c
c     # Limiters are applied to the fwave's, and the only significant
c     # modification of this code is in the "do 110" loop, for the 
c     # second order corrections.
c     ----------------------------------------------------------------------
c
c     method(1) = 1   ==>  Godunov method
c     method(1) = 2   ==>  Slope limiter method
c     mthlim(p)  controls what limiter is used in the pth family
c
c
c     amdq, apdq, fwave, s, and f are used locally:
c
c     amdq(1-mbc:maxmx+mbc, meqn) = left-going flux-differences
c     apdq(1-mbc:maxmx+mbc, meqn) = right-going flux-differences
c        e.g. amdq(i,m) = m'th component of A^- \Delta q from i'th Riemann
c                         problem (between cells i-1 and i).
c
c     fwave(1-mbc:maxmx+mbc, meqn, mwaves) = waves from solution of
c                                           Riemann problems,
c            fwave(i,m,mw) = mth component of jump in f across
c                           wave in family mw in Riemann problem between
c                           states i-1 and i.
c
c     s(1-mbc:maxmx+mbc, mwaves) = wave speeds,
c            s(i,mw) = speed of wave in family mw in Riemann problem between
c                      states i-1 and i.
c
c     f(1-mbc:maxmx+mbc, meqn) = correction fluxes for second order method
c            f(i,m) = mth component of flux at left edge of ith cell 
c     --------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension    q(1-mbc:maxmx+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, *)
      dimension    f(1-mbc:maxmx+mbc, meqn)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension fwave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension amdq(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)
      dimension dtdx(1-mbc:maxmx+mbc)
      dimension method(7),mthlim(mwaves)
      logical limit
c
c     # check if any limiters are used:
      limit = .false.
      do 5 mw=1,mwaves
	 if (mthlim(mw) .gt. 0) limit = .true.
   5     continue
c
      mcapa = method(6)
      do 10 i=1-mbc,mx+mbc
	 if (mcapa.gt.0) then
	     if (aux(i,mcapa) .le. 0.d0) then
		write(6,*) 'Error -- capa must be positive'
		stop
		endif
             dtdx(i) = dt / (dx*aux(i,mcapa))
	    else
             dtdx(i) = dt/dx
	    endif
   10	 continue      
c
c
c
c     # solve Riemann problem at each interface 
c     -----------------------------------------
c
      call rp1(maxmx,meqn,mwaves,mbc,mx,q,q,aux,aux,fwave,s,amdq,apdq)
c      stop
c
c   LOOP FOR DEBUG
c      
c         do i = 1-mbc,mx+mbc
c            print*, q(i,1), q(i,2), q(i,3), q(i,4), q(i,5), q(i,6)
c            print*, q(i,4), q(i,6), q(i,5), q(i,2), q(i,7), q(i,8)
c         enddo
c
c     # Modify q for Godunov update:
c     # Note this may not correspond to a conservative flux-differencing
c     # for equations not in conservation form.  It is conservative if
c     # amdq + apdq = f(q(i)) - f(q(i-1)).
c
      do i=1,mx+1
         do m=1,meqn
            q(i,m)   = q(i,m)   - dtdx(i)*apdq(i,m)
            q(i-1,m) = q(i-1,m) - dtdx(i-1)*amdq(i,m)
         enddo
      enddo
c
c   LOOP FOR DEBUG
c       
c         do i = 1-mbc,mx+1
c             print*, apdq(i,1), apdq(i,2), apdq(i,3), apdq(i,4),
c     &               apdq(i,5), apdq(i,6)
c             print*, apdq(i,4), apdq(i,6), apdq(i,5), apdq(i,2),
c     &               apdq(i,7), apdq(i,8)
c            print*, q(i,1), q(i,2), q(i,3), q(i,4), q(i,5), q(i,6)
c            print*, q(i,4), q(i,6), q(i,5), q(i,2), q(i,7), q(i,8)
c         enddo
         
c
c     # compute maximum wave speed:
      cfl = 0.d0
      do 50 mw=1,mwaves
	 do 45 i=1,mx+1
c          # if s>0 use dtdx(i) to compute CFL,
c          # if s<0 use dtdx(i-1) to compute CFL:
	   cfl = dmax1(cfl, dtdx(i)*s(i,mw), -dtdx(i-1)*s(i,mw))
   45      continue
   50    continue
c
      if (method(2) .eq. 1) go to 900
c
c     # compute correction fluxes for second order q_{xx} terms:
c     ----------------------------------------------------------
c     
      do 100 m = 1, meqn
            do 100 i = 1-mbc, mx+mbc
               f(i,m) = 0.d0
  100          continue
c
c     # apply limiter to waves:
      if (limit) then
         call limiter(maxmx,meqn,mwaves,mbc,mx,fwave,s,mthlim)
         endif

c
      do 120 i=1,mx+1
	 do 120 m=1,meqn
	    do 110 mw=1,mwaves
	       dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
	       f(i,m) = f(i,m) + 0.5d0 * dsign(1.d0,s(i,mw))
     &		   * (1.d0 - dabs(s(i,mw))*dtdxave) * fwave(i,m,mw)
c

  110          continue
  120       continue
c
c
  140 continue
c
c     # update q by differencing correction fluxes 
c     ============================================
c
c     # (Note:  Godunov update has already been performed above)
c
      do 150 m=1,meqn
	 do 150 i=1,mx
	    q(i,m) = q(i,m) - dtdx(i) * (f(i+1,m) - f(i,m))
  150       continue
c
  900 continue
      return
      end
