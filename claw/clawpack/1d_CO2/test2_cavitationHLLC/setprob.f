      subroutine setprob
c
c
      USE properties
c
      implicit double precision (a-h,o-z)
      common /comic/ rhol,rhoul,el,rhor,rhour,er,Tl,Tr
      common /e_cst/ e_const 
c
c     # Set gamma and gamma1 = gamma-1 for Euler equations
c     # Passed to the Riemann solver rp1.f in a common block

c     # set values for use in qinit: shock tube initial data
c
      e_const = 0.5d6 
      open(unit=7,file='setprob.data',status='old',form='formatted')


       read(7,*) rhol,ul,Tl
       read(7,*) rhor,ur,Tr
c
c      # data in left state:
       rhoul = rhol*ul
       v = 1.0e0 / rhol
       
       CALL inter_energy(Tl,v,e)
       el = e+e_const
       el = rhol*(el + 0.5d0*ul*ul)
c
c      # data in right state:
       rhour = rhor*ur
       v = 1.0e0 / rhor
       CALL inter_energy(Tr,v,e)
       er = e+e_const
       er = rhor*(er + 0.5d0*ur*ur)
c
      return
      end
