c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dxc,dyc,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays for advection on a curvilinear grid
c
c     # on input, (xc(i),yc(j)) gives uniformly spaced computational grid.
c     # on output, 
c     #   aux(i,j,1) is nx at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is ny at "left" boundary of grid point (i,j)
c     #   aux(i,j,3) is gamma ratio at "left" boundary of grid point (i,j)
c     #   aux(i,j,4) is nx at "bottom" boundary of grid point (i,j)
c     #   aux(i,j,5) is ny at "bottom" boundary of grid point (i,j)
c     #   aux(i,j,6) is gamma ratio at "bottom" boundary of grid point (i,j)
c     #   aux(i,j,7) = kappa  is ratio of cell area to dxc*dyc
c     
      implicit double precision (a-h,o-z)
      double precision, intent(out), 
     &     dimension (1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux) :: aux
c      dimension  aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      dimension xccorn(5),yccorn(5),xpcorn(5),ypcorn(5)
c
c
c      print*, 'setaux.f starts'
      do 20 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
c
c           # computational points (xc,yc) are mapped to physical
c           # coordinates (xp,yp) by mapc2p:
c
c           # lower left corner
            xccorn(1) = xlower + (i-1)*dxc
            yccorn(1) = ylower + (j-1)*dyc
            call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1))
!            print*, xccorn(1),yccorn(1),xpcorn(1),ypcorn(1)

c           # upper left corner:
            xccorn(2) = xccorn(1)
            yccorn(2) = yccorn(1) + dyc
            call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2))
c
c           # upper right corner:
            xccorn(3) = xccorn(1) + dxc
            yccorn(3) = yccorn(1) + dyc
            call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3))
c
c           # lower right corner:
            xccorn(4) = xccorn(1) + dxc
            yccorn(4) = yccorn(1)
            call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4))
!            print*, xccorn(4),yccorn(4),xpcorn(4),ypcorn(4)
c           
            xpcorn(5) = xpcorn(1)
            ypcorn(5) = ypcorn(1)
c
            xccorn(5) = xccorn(1)
            yccorn(5) = yccorn(1)
c
c           # normal vectors
c
            dy1 = ypcorn(2) - ypcorn(1)
            dx1 = xpcorn(2) - xpcorn(1)
            p_norm1 = sqrt(dy1*dy1 + dx1*dx1)
            c_norm1 = sqrt( (yccorn(2) - yccorn(1))**(2.0) +
     &                      (xccorn(2) - xccorn(1))**(2.0) )
c
c            print*,'nx', dy1/p_norm
c            print*,'ny', dx1/p_norm
            aux(i,j,1) =  dy1/p_norm1
            aux(i,j,2) = -dx1/p_norm1
            aux(i,j,3) = p_norm1/dyc 
c            print*,i,j,aux(i,j,1),aux(i,j,2),aux(i,j,3)
c            print*,'gamma1',i,j,aux(i,j,3)
c
            dy4 = ypcorn(5) - ypcorn(4)
            dx4 = xpcorn(5) - xpcorn(4)
            p_norm4 = sqrt(dy4*dy4 + dx4*dx4)
            c_norm4 = sqrt( (yccorn(1) - yccorn(4))**(2.0) + 
     &                      (xccorn(1) - xccorn(4))**(2.0) )
c
            aux(i,j,4) =   dy4/p_norm4     
            aux(i,j,5) =  -dx4/p_norm4    
            aux(i,j,6) = p_norm4/dxc  
c            print*,'gamma4',aux(i,j,6)

c
c
c           # compute area of physical cell from four corners:

            area = 0.d0  
            do ic=1,4
            area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) *
     &               (xpcorn(ic+1)-xpcorn(ic)) 
            enddo
	    
            aux(i,j,7) = area / (dxc*dyc)
c
   20  continue
       return

       end
