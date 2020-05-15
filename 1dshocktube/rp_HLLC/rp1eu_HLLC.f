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
c
c      USE def_variables
c      USE def_constants
c      USE properties
c      USE interp_functions
c      USE non_linear_solvers
      USE Interp_table
c
c
      implicit double precision (a-h,o-z)
c
c      
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
      dimension qml(3), qmr(3)
      real(8)  press, sound, vr, vl, press_guess, 
     &        xr, xl, a_out, res,Tl, Tr, ar, al
      real(8)  press_r, press_l, cl, cr
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
         vl = 1.0d0 / rl
c
         rr = qr(i,  1)
         vr = 1.0d0 / rr
c
         ul = ql(i-1,2)/rl
         ur = qr(i,  2)/rr
         el = ql(i-1,3)/rl - ul*ul/2d0 - 5.0d6
         er = qr(i,  3)/rr - ur*ur/2d0 - 5.0d6
c
c        values from TABLE CO2
         CALL CO2BLLT_EQUI(press, Tl, sound, xl, al, res,
     &                 el,vl,p_guess)
c         
         press_l = press
         cl = sound
         
             
         CALL CO2BLLT_EQUI(press, Tr, sound, xr, ar, res,
     &                 er,vr,p_guess)
         press_r = press
         cr = sound
c 
c        speed of sound  for signal velocity calculation
c 
         sl     = min(ul-cl, ur-cr)
         sr     = max(ul+cl, ur+cr)
         sm     = (press_r-press_l + rl*ul*(sl- ul) - rr*ur*(sr- ur))
     &                  / (rl*   (sl- ul) - rr*   (sr- ur))
         s(i,1) = sl
         s(i,2) = sm
         s(i,3) = sr
c
c        intermediate states
c
         qml(1) = 1d0
         qml(2) = sm
         qml(3) = (ql(i-1,3)/rl + (sm - ul)*(sm + press_l/
     &                        (rl*(sl - ul))))
         qml(1:3) = qml(1:3) * rl*(sl - ul)/(sl - sm) 
c
         qmr(1) = 1d0
         qmr(2) = sm
         qmr(3) = (qr(i,  3)/rr + (sm - ur)*(sm + press_r/
     &                        (rr*(sr - ur))))
         qmr(1:3) = qmr(1:3) * rr*(sr - ur)/(sr - sm) 
c
c        # Compute the waves.
c
         wave(i,1:3,1) = qml(1:3)  - ql(i-1,1:3)
         wave(i,1:3,2) = qmr(1:3)  - qml(1:3) 
         wave(i,1:3,3) = qr(i,1:3) - qmr(1:3) 
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
