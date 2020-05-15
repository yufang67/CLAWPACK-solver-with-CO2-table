c
c
c =========================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q, not ghost points.
c     # shocktube with states set in setprob.f
c
c
      USE Grid
      USE properties, ONLY: inter_energy
c      USE peng_robinson, ONLY: interenergy_pr
c      USE stiffened
      USE var_const, ONLY: e_const, guessP,eguess_out,
     &                     vguess_in,Tguess


      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
      common /comic/ rhol,rhoul,el,rhor,rhour,er,Tl,Tr,rhovl,rhovr
      common /flags/ i_flag
c
c
      CALL cpu_time(time_1)
      CALL MAKE_GRID()
      CALL cpu_time(time_2)

c
c     
         print*, "Span-Wagner EoS"
       ul = rhoul / rhol
       velol = rhovl / rhol
c
       CALL inter_energy(Tl,1.0/rhol,eint)
       el = eint + e_const
       el = rhol*(el + 0.5d0*(ul*ul+velol*velol))
c       el = rhol*(e + 0.5d0*(ul*ul+velol*velol))
c       el = el+e_const
c
c         print*,'left state', eint, e_const
       ur = rhour / rhor
       velor = rhovr / rhor
c
       CALL inter_energy(Tr,1.0/rhor,eint)
       er = eint + e_const
       er = rhor*(er + 0.5d0*(ur*ur+velor*velor))
c       er = rhor*(e + 0.5d0*(ur*ur+velor*velor))
c       er = er+e_const
c         print*,'right state', eint, e_const

c
c  to adapt to 2D 
       sloc = 50.0
c       sloc = 200.0
c
      do 150 j=1-mbc,my+mbc
c          guessp_BC(j) = 1.0d6
c          vguess_in(j)  = 1.0/500.0
c          eguess_out(j) = -108.58d3    
          do 150 i=1-mbc,mx+mbc
             xcell = xlower + (i-0.5d0)*dx
             if (xcell .lt. sloc) then
                 q(i,j,1) = rhol
                 q(i,j,2) = rhoul
                 q(i,j,3) = rhovl
                 q(i,j,4) = el
c                 guessP_rpn(i,j) = 0.19d6
                 guessP(i,j) = 9d6
c                 Tguess(i,j) = 400
c                 c_BC(i,j)   = 194.0
              else
                 q(i,j,1) = rhor
                 q(i,j,2) = rhour
                 q(i,j,3) = rhovr
                 q(i,j,4) = er
c                 guessP_rpn(i,j) = 0.15d6
                 guessP(i,j) = 7d6
c                 Tguess(i,j) = 400
c                 c_BC(i,j)   = 239.0
              endif
c        print*, q(i,j,1), q(i,j,2), q(i,j,3), q(i,j,4)
  150        continue
      return
      end
