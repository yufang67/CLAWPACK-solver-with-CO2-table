c
c
c
c =========================================================
      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &           wave,s,amdq,apdq)
c =========================================================
c
c     # solve Riemann problems for the 1D Euler equations using Roe's 
c     # approximate Riemann solver.  
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit double precision (a-h,o-z)
      dimension   ql(1-mbc:maxmx+mbc, meqn)
      dimension   qr(1-mbc:maxmx+mbc, meqn)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension amdq(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)
c
c     # local storage
c     ---------------
      parameter (max2 = 2002)  !# assumes at most 2000 grid points with mbc=2
      dimension Fl(3), Fr(3), qm(3)
      dimension u(-1:max2),enth(-1:max2),a(-1:max2)
      common /param/  gamma,gamma1
c
c
c     # Compute Roe-averaged quantities and fluxes
c
c
      do 20 i=2-mbc,mx+mbc
c
c        left and right properties
c
         rl = ql(i-1,1)
         rr = qr(i,  1)
         ul = ql(i-1,2)/rl
         ur = qr(i,  2)/rr
         el = ql(i-1,3)/rl - ul*ul/2d0
         er = qr(i,  3)/rr - ur*ur/2d0
         pl = gamma1*rl*el                  ! p = (gamma-1) * rho e 
         pr = gamma1*rr*er
         cl = sqrt(gamma * pl / rl)         ! c = (gamma * p / rho)^(0.5) 
         cr = sqrt(gamma * pr / rr)
c 
c        Roe-averaged quantities:
c
         rhsqrtl = sqrt(rl)
         rhsqrtr = sqrt(rr)
         enth(i) = (((ql(i-1,3)+pl)/rhsqrtl
     &             + (qr(i,  3)+pr)/rhsqrtr)) / (rhsqrtl + rhsqrtr)
c        Roe-average particle speed
         u(i) = (rhsqrtl*ul + rhsqrtr*ur) / (rhsqrtl + rhsqrtr)     
c        Roe-average sound speed
         a(i) = gamma1*(enth(i) - .5d0*u(i)**2)
         a(i) = sqrt(a(i))
c
c        speed of sound  for signal velocity calculation
c 
         s(i,1) = min(ul-cl, u(i)-a(i), 0d0)
         s(i,2) = max(ur+cr, u(i)+a(i), 0d0)
c
c        FLUXES
c
         Fl(1) = rl * ul
         Fl(2) = rl * ul*ul + pl
         Fl(3) = (ql(i-1,3) + pl)* ul
         Fr(1) = rr * ur
         Fr(2) = rr * ur*ur + pr
         Fr(3) = (qr(i,  3) + pr)* ur
c
c        intermediate state
c
         qm(1:3) = ( s(i,2)*qr(i,1:3) - s(i,1)*ql(i-1,1:3)
     &             - Fr(1:3) + Fl(1:3)) / (s(i,2) - s(i,1))
c
c        # Compute the waves.
c
         wave(i,1:3,1) = qm(1:3)   - ql(i-1,1:3)
         wave(i,1:3,2) = qr(i,1:3) - qm(1:3) 
   20    continue
c
c     # compute Godunov flux f0:
c     --------------------------
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do 100 m=1,3
         do 100 i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do 90 mw=1,mwaves
               if (s(i,mw) .lt. 0.d0) then
                   amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
                 else
                   apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                 endif
   90          continue
  100       continue
      return
      end
