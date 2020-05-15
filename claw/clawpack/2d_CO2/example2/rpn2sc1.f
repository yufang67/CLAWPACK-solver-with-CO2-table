
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the sample scalar equation
c     #  u_t + (u^2)_x + (u^4)_y = 0
c     
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1 
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q, 
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #               
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c     # maux=0 and aux arrays are unused in this example.
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension  apdq(1-mbc:maxm+mbc, meqn)
      dimension  amdq(1-mbc:maxm+mbc, meqn)
c
      logical efix
      parameter (maxm2 = 502)
      common /comrp/ stran(-1:maxm2)
c
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
c
      if (ixy.eq.1) then
c         # Solve Riemann problem in x-direction for u_t + (u^2)_x = 0:
c           -----------------------------------------------------------
c
          do 10 i = 2-mbc, mx+mbc
c            # wave is jump in q, speed comes from R-H condition:
             wave(i,1,1) = ql(i,1) - qr(i-1,1)
             s(i,1) = qr(i-1,1) + ql(i,1) 
c
c            # compute left-going and right-going flux differences:
c
             df = ql(i,1)**2 - qr(i-1,1)**2
             if (s(i,1) .lt. 0.d0) then
                 amdq(i,1) = df
                 apdq(i,1) = 0.d0
               else
                 amdq(i,1) = 0.d0
                 apdq(i,1) = df
               endif
c
c            # check for sonic point, in which case f(q^0) = 0:
             if (efix .and. qr(i-1,1).lt.0.d0 
     &                .and. ql(i,1).gt.0.d0) then
                 amdq(i,1) = -qr(i-1,1)**2
                 apdq(i,1) = ql(i,1)**2
                 endif
c
c            # set velocity in transverse direction, passed to rpt2:
             stran(i) = qr(i-1,1)**3 + ql(i,1)*qr(i-1,1)**2
     &               + ql(i,1)**2*qr(i-1,1) + ql(i,1)**3
   10        continue
c
c
          else  !# (ixy = 2)
c         # Solve Riemann problem in y-direction for u_t + (u^4)_y = 0:
c           -----------------------------------------------------------
c
          do 20 i = 2-mbc, mx+mbc
c            # wave is jump in q, speed comes from R-H condition:
             wave(i,1,1) = ql(i,1) - qr(i-1,1)
             s(i,1) = qr(i-1,1)**3 + ql(i,1)*qr(i-1,1)**2
     &               + ql(i,1)**2*qr(i-1,1) + ql(i,1)**3
c
c            # compute left-going and right-going flux differences:
c
             df = ql(i,1)**4 - qr(i-1,1)**4
             if (s(i,1) .lt. 0.d0) then
                 amdq(i,1) = df
                 apdq(i,1) = 0.d0
               else
                 amdq(i,1) = 0.d0
                 apdq(i,1) = df
               endif
c
c            # check for sonic point, in which case f(q^0) = 0:
             if (efix .and. qr(i-1,1).lt.0.d0 
     &                .and. ql(i,1).gt.0.d0) then
                 amdq(i,1) = -qr(i-1,1)**4
                 apdq(i,1) = ql(i,1)**4
                 endif
c
c            # set velocity in transverse direction, passed to rpt2:
             stran(i) = qr(i-1,1) + ql(i,1) 
   20        continue
         endif
      return
      end
