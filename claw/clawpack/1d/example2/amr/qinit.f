c
c
c =========================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c
c     # 1d data written to initialize a 2d grid for use with amrclaw
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,*)
      common /cqinit/ beta
c
c
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         q(i,1,1) = dexp(-50 * (xcell-0.3d0)**2)  
     &            * dcos(20.d0*(xcell-0.3d0))
         if (xcell .gt. 0.30d0) q(i,1,1) = 0.d0
         q(i,1,2) = 0.d0
  150    continue
c
      return
      end
