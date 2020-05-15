c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # shocktube with states set in setprob.f
c
c
      USE Grid
      USE properties
      USE peng_robinson, ONLY: interenergy_pr
      USE stiffened


      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /comic/ rhol,rhoul,el,rhor,rhour,er,Tl,Tr
      common /flags/ i_flag
      common /e_cst/ e_const
c
c     
c     print*, i_flag
      CALL cpu_time(time_1)
      CALL MAKE_GRID()
      CALL cpu_time(time_2)
      Print*, 'Make table grid takes', time_2-time_1
c 
c
      IF (i_flag==1) THEN
       print*, "1 Span-Wagner LOOK-UP TABLE"
       v = 1.0e0 / rhol
       ul = rhoul / rhol
c
       CALL inter_energy(Tl,v,e)
       el = e+e_const
       el = rhol*(el + 0.5d0*ul*ul)
c
       v = 1.0e0 / rhor
       ur = rhour / rhor
c
       CALL inter_energy(Tr,v,e)
       er = e+e_const
       er = rhor*(er + 0.5d0*ur*ur)
c
c     Create TABLE CO2     
c      CALL MAKE_GRID()
c
c
c
      ELSEIF (i_flag==2)THEN
      print*, "2 Peng-robinson EOS"
c     
      v = 1.0e0 / rhol
      ul = rhoul / rhol
c     
      CALL interenergy_pr(Tl,v,el)
c      print*, Tl,v,el
      el = rhol*(el + 0.5d0*ul*ul)
c
      v = 1.0e0 / rhor
      ur = rhour / rhor
c     
      CALL interenergy_pr(Tr,v,er)
      er = rhor*(er + 0.5d0*ur*ur)
c
c
c
      ELSEIF (i_flag==3)THEN
      print*, "3 Stiffened gas EOS"
c
      ul = rhoul / rhol
      ur = rhour / rhor
      el = el+e_const
      er = er+e_const
      el = rhol*(el + 0.5d0*ul*ul)
      er = rhor*(er + 0.5d0*ur*ur)
c
c
c
      ELSEIF (i_flag==4)THEN
      print*, "4 Span-wagner EOS"
c 
       v = 1.0e0 / rhol
       ul = rhoul / rhol
c
       CALL inter_energy(Tl,v,e)
       el = e+e_const
       el = rhol*(el + 0.5d0*ul*ul)
c
       v = 1.0e0 / rhor
       ur = rhour / rhor
c
       CALL inter_energy(Tr,v,e)
       er = e+e_const
       er = rhor*(er + 0.5d0*ur*ur)
c      
c
c
      ENDIF
c
c
c 
       sloc = 0.5d0
c
       do 150 i=1,mx
             xcell = xlower + (i-0.5d0)*dx
             if (xcell .lt. sloc) then
                 q(i,1) = rhol
                 q(i,2) = rhoul
                 q(i,3) = el
                else
                 q(i,1) = rhor
                 q(i,2) = rhour
                 q(i,3) = er
                endif
  150        continue
      
      return
      end
