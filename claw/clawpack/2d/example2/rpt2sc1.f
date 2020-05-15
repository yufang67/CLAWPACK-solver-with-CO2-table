c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the scalar equation
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # For the scalar equation, this simply amounts to computing the
c     # transverse wave speed sb from the transverse Riemann problem and
c     # setting bmasdq or bpasdq to sb*asdq.
c
      dimension     ql(1-mbc:maxm+mbc, meqn)
      dimension     qr(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
      parameter (maxm2 = 502)
      common /comrp/ stran(-1:maxm2)
c
c     # transverse wave speeds have been computed in rpn2
c     # maux=0 and aux arrays are unused in this example.
c
      do 10 i = 2-mbc, mx+mbc
          bmasdq(i,1) = dmin1(stran(i), 0.d0) * asdq(i,1)
          bpasdq(i,1) = dmax1(stran(i), 0.d0) * asdq(i,1)
   10     continue
c
      return
      end
