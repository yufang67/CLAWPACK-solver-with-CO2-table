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
c      USE properties
c      USE interp_functions
c      USE non_linear_solvers
      USE Interp_table
      USE peng_robinson
      USE solver_eos
      USE location
      USE stiffened
      USE properties
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
      common /e_cst/ e_const,Tguess
      common /flags/ i_flag
c
c     # local storage
c     ---------------
      parameter (max2 = 2002)  !# assumes at most 2000 grid points with mbc=2
      dimension qml(3), qmr(3)
      integer exitflag
      real(8)  press, sound, vr, vl, 
     &        xr, xl, a_out, res,Tl, Tr, ar, al
      real(8)  press_r, press_l, cl, cr, guess_t
c
c
c     # Compute Roe-averaged quantities and fluxes
c
c
c     print*, Tguess
      do 20 i=2-mbc,mx+mbc
c
c left and right properties
c
         rl = ql(i-1,1)
         vl = 1.0d0 / rl
c
         rr = qr(i,  1)
         vr = 1.0d0 / rr
c
         ul = ql(i-1,2)/rl
         ur = qr(i,  2)/rr
         el = ql(i-1,3)/rl - ul*ul/2d0
         er = qr(i,  3)/rr - ur*ur/2d0
c
c  EOS:  Using er,el, rhor, rhol to find press_r and press_l,c_r,c_l
      IF (i_flag==1)THEN
c
c  values from TABLE CO2
c
         pguess = 2.0e6
         el = el - e_const
         CALL CO2BLLT_EQUI(press, Tl, sound, xl, al, res1,
     &                 el,vl,pguess)
c        
         press_l = press
         cl = sound
         pguess = press
         
         er = er - e_const    
         CALL CO2BLLT_EQUI(press, Tr, sound, xr, ar, res2,
     &                 er,vr,pguess)
         press_r = press
         cr = sound
         pguess = press
       
!        IF (i==500 .OR. i==501) THEN
!     print*,i, res1, res2
!        ENDIF
       
      ELSE IF (i_flag==2) THEN
         energy = el - 241.3d3
         CALL phaseloca(xl, al, i_flag_phase, energy,vl)
         IF (i_flag_phase==1 .OR. i_flag_phase==5) THEN
            PRINT*, 'liquid presents in the left wave in rp1.f at',i
c            STOP
         ENDIF
c
         IF (i.EQ.0) guess_t=300.0+1.0
         CALL eos_1d(2, Tl, out_2, resnorm, Niter,
     &            exitflag, el, guess_t, vl, out3)
         guess_t = Tl+1.0
         CALL pressure_pr(Tl,vl,press)
         CALL soundspeed_pr(Tl,vl,sound)
         press_l = press
         cl = sound
c     print*, "i node left",i, Tl,press_l,cl
         energy = er - 241.3d3
         CALL phaseloca(xr, ar, i_flag_phase, energy,vr)
         IF (i_flag_phase==1 .OR. i_flag_phase==5) THEN
            PRINT*, 'liquid presents in the right wave in rp1.f at',i
c            STOP
         ENDIF
c
         CALL eos_1d(2, Tr, out_2, resnorm, Niter,
     &            exitflag, er, guess_t, vr, out3)
         guess_t = Tr-1.0
         CALL pressure_pr(Tr,vr,press)
         CALL soundspeed_pr(Tr,vr,sound)
         press_r = press
         cr = sound
!     print*, 'temperature', Tl, Tr, guess_t
c     print*, "i node right",i, Tr,press_r,cr
      ELSE IF (i_flag==3) THEN
         energy = el - e_const
         CALL phaseloca(xl, al, i_flag_phase, energy,vl)
c
c      print*, 'flag_phase=', i_flag_phase,i
         CALL pressure_st(i_flag_phase,vl,energy, press_l)
         CALL sound_speed_st(i_flag_phase,vl,energy,cl)
c
         energy = er - e_const
         CALL phaseloca(xr, ar, i_flag_phase, energy,vr)
c
c      print*, 'flag_phase=', i_flag_phase,i
         CALL pressure_st(i_flag_phase,vr,energy, press_r)
         CALL sound_speed_st(i_flag_phase,vr,energy,cr)
c
c
      ELSE IF (i_flag==4) THEN
         energy = el - e_const
         CALL phaseloca(xl, al, i_flag_phase, energy,vl)
c     
        IF (i_flag_phase==5) THEN
           IF (i.EQ.1) pguess_out = 6.0e6
           CALL CO2BLLT_EQUI(press_l, Tl, cl, xl, al, res1,
     &                     energy, vl, pguess_out)
           pguess_out = press_l
           i_flag_phase = res1
        ELSE
c         IF (i.EQ.1) Tguess=Tl
           CALL eos_1d(4, Tl, out_2, resnorm, Niter,
     &            exitflag, energy, Tguess, vl, out3)
           Tguess = Tl
           CALL pressure(Tl,vl,press_l)
           CALL sound_speed(Tl,vl,cl)
        ENDIF
c
        energy = er - e_const
         CALL phaseloca(xr, ar, i_flag_phase, energy,vr)
c     
        IF (i_flag_phase==5) THEN
           IF (i.EQ.1) pguess_out = 6.0e6
           CALL CO2BLLT_EQUI(press_r, Tr, cr, xr, ar, res2,
     &                     energy, vr, pguess_out)
           pguess_out = press_r
           i_flag_phase = res2
        ELSE
c         IF (i.EQ.1) Tguess=Tl
           CALL eos_1d(4, Tr, out_2, resnorm, Niter,
     &            exitflag, energy, Tguess, vr, out3)
           Tguess = Tr
           CALL pressure(Tr,vr,press_r)
           CALL sound_speed(Tr,vr,cr)
        ENDIF
      
      ENDIF
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
