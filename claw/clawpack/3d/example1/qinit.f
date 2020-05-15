
c
c
c
c     =====================================================
      subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &     xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
c
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &     1-mbc:maxmz+mbc, meqn)

c
c     # set concentration profile
c     ---------------------------
      do k=1,mz
         zcell = zlower + (k-0.5d0)*dz
         do j=1,my
            ycell = ylower + (j-0.5d0)*dy
            do i=1,mx
               xcell = xlower + (i-0.5d0)*dx
               if (xcell .gt. 0.1d0 .and. xcell .lt. 0.5d0 .and.
     &               ycell .gt. 0.1d0 .and. ycell .lt. 0.5d0 .and.
     &               zcell .gt. 0.1d0 .and. zcell .lt. 0.5d0) then
                  q(i,j,k,1) = 1.d0
               else
                  q(i,j,k,1) = 0.d0
               endif
            end do
         end do
      enddo

      return
      end
