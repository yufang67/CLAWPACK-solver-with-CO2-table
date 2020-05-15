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
c     # parameters of the nozzle
      convL = 40
      divL = 60
      angle_c = 0.09348
      angle_d = 0.05794
      x_t = 40
      y_t = 2.5/2
c      print*, convL,divL, angle_c,angle_d,x_t,y_t
c     stop
c     ##############################
      x_f = convL+divL
      y_s = convL*tan(angle_c) + y_t
      y_f = divL*tan(angle_d) + y_t
c      
      AA_up   = (y_t - y_s) / x_t  
      AA_down = (y_s - y_t) / x_t
c
      BB_up   = (y_f - y_t) / (x_f - x_t)
      BB_down = (y_t - y_f) / (x_f - x_t) 
c     ##################################
c
c      xc = xc*0.0835
      IF (xc .LE. x_t) THEN 
         IF (yc .GE. y_s) THEN
            yp = (xc *AA_up + y_s)*(yc-y_s)/y_s
            xp = xc
         ELSE
            yp = (xc *AA_down - y_s)*(y_s-yc)/y_s  
            xp = xc
         ENDIF
      ELSE
         IF (yc .GE. y_s) THEN
            yp = (xc *BB_up + y_f - BB_up*x_f)*(yc-y_s)/y_s
            xp = xc
         ELSE
            yp = (xc *BB_down -  y_f - x_f*BB_down)*(y_s-yc)/y_s
            xp = xc
         ENDIF
      ENDIF
      
      return
      end
