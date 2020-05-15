      subroutine setprob
c
c
      USE properties
c
      implicit double precision (a-h,o-z)
c      integer flag_eos
      common /comic/ rhol,rhoul,el,rhor,rhour,er,Tl,Tr
      common /e_cst/ e_const,Tguess
      common /flags/ i_flag
c
c     # Passed to the Riemann solver rp1.f in a common block

c     # set values for use in qinit: shock tube initial data
c
      e_const = 0.5d6 
      open(unit=7,file='setprob.data',status='old',form='formatted')

       read(7,*) i_flag
       IF (i_flag==3) THEN
       read(7,*) rhol,ul,Tl,el
       read(7,*) rhor,ur,Tr,er
       ELSE
       read(7,*) rhol,ul,Tl
       read(7,*) rhor,ur,Tr
       ENDIF
c
       Tguess = Tl+3.0
c      # data in left state:
c       print*, i_flag       

c      IF (flag_eos==1) THEN
c       print*, "EOS is CO2 LOOK-UP TABLE"
       rhoul = rhol*ul
c       v = 1.0e0 / rhol
       
c       CALL inter_energy(Tl,v,e)
c       el = e+e_const
c       el = rhol*(el + 0.5d0*ul*ul)
c
c      # data in right state:
       rhour = rhor*ur
c       v = 1.0e0 / rhor
c       CALL inter_energy(Tr,v,e)
c       er = e+e_const
c       er = rhor*(er + 0.5d0*ur*ur)
c      ENDIF
c
      return
      end
