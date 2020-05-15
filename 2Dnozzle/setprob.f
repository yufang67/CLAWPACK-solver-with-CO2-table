      subroutine setprob
c
c
c      USE properties
c      USE var_const, ONLY: e_const,i_flag
c
      implicit double precision (a-h,o-z)
      common /comic/ rhol,rhoul,el,rhor,rhour,er,Tl,Tr,rhovl,rhovr
      common /flags/ i_flag 
c
c     # Passed to the Riemann solver rp1.f in a common block

c     # set values for use in qinit: shock tube initial data
c
c      e_const = 0.5d6 
      open(unit=7,file='setprob.data',status='old',form='formatted')

       read(7,*) i_flag
       IF (i_flag==3) THEN
       read(7,*) rhol,ul,velol,Tl,el
       read(7,*) rhor,ur,velor,Tr,er
       ELSE
       read(7,*) rhol,ul,velol,Tl
       read(7,*) rhor,ur,velor,Tr
       ENDIF
c
       rhoul = rhol*ul
       rhovl = rhol*velol
c
c      # data in right state:
       rhour = rhor*ur
       rhovr = rhor*velor
c       print*, 'read finished'
      return
      end
