c
c     =====================================================
      subroutine mapc2p(xc,yc,xp,yp)
c     =====================================================
c
c     # on input,  (xc,yc) is a computational grid point
c     # on output, (xp,yp) is corresponding point in physical space
c
      implicit double precision (a-h,o-z)
c
      pi = 3.1415926
c     # parameters of the nozzle09_noz1
      convL = 27.35d-3
      divL  = 56.15d-3
      thL   = 0.0
      sp1L  = 2.0d-3
      sp2L  = 2.0d-3
      ttan  = 1.32645d-3 
c     #
      thW   = 0.24d-3
      axi   = 10.0
      axi_py= 5.0d-3
      x1    = sp1L
      x2    = sp1L+convL
      x3    = sp1L+convL+thL
      x4    = sp1L+convL+thL+divL
      x5    = sp1L+convL+thL+divL+sp2L
c      print*, convL,divL, angle_c,angle_d,x_t,y_t
c     stop
c     ##############################
      AA   =  (axi_py-thW/2.0)/convL 
c      BB   =  (axi_py-thW/2.0)/divL
      BB   =  ttan
c     ##################################
c
      xp = xc /100.0*(sp1L+convL+thL+divL+sp2L) 
c
c     #################################
c
      IF (xp .LE. x1) THEN
c## Sponge
         yp  =  (yc)/(axi)*axi_py
      ELSEIF ((xp .GE. x1) .AND. (xp .LE. x2)) THEN
c## Convergence
         IF (yc .GE. axi) THEN
            yp = ((x2-xp)*AA+thW/2.0)*(yc-axi)/axi + axi_py
         ELSE
            yp = axi_py - ((x2-xp)*AA+thW/2.0)*(axi-yc)/axi
         ENDIF
c
      ELSEIF ((xp .GT. x2 ) .AND. (xp .LE. x3)) THEN
c##throat
         IF (yc .GE. axi) THEN
            yp = (yc-axi)/axi*thW/2.0  + axi_py
         ELSE 
            yp = axi_py - (axi-yc)/axi*thW/2.0
         ENDIF
c
      ELSEIF ((xp .GT. x3) .AND. (xp .LE. x4)) THEN
c##divergence
         IF (yc .GE. axi) THEN
            yp = ((divL-x4+xp)*BB+thW/2.0)*(yc-axi)/axi + axi_py
         ELSE
            yp = axi_py - ((divL-x4+xp)*BB+thW/2.0)*(axi-yc)/axi
         ENDIF
c##sponge
      ELSEIF ((xp .GE. x4)) THEN
         IF (yc .GE. axi) THEN
            yp  =  (yc-axi)/axi*(divL*BB+thW/2.0) + axi_py
         ELSE
            yp  =  axi_py - (axi-yc)/axi*(divL*BB+thW/2.0)
         ENDIF
      ENDIF
c      
      return
      end
