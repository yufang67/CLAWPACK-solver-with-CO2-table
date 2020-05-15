      subroutine setprob
      implicit double precision (a-h,o-z)
      common /param/ gamma, gamma1
      common /comic/ rhol,rhoul,el,rhor,rhour,er
c
c     # Set gamma and gamma1 = gamma-1 for Euler equations
c     # Passed to the Riemann solver rp1.f in a common block

c     # set values for use in qinit: shock tube initial data
c
      open(unit=7,file='setprob.data',status='old',form='formatted')


      read(7,*) gamma
      gamma1 = gamma - 1.d0

       read(7,*) rhol,ul,pl
       read(7,*) rhor,ur,pr
c
c      # data in left state:
       rhoul = rhol*ul
       el = pl/gamma1 + 0.5d0*rhol*ul**2
c
c      # data in right state:
       rhour = rhor*ur
       er = pr/gamma1 + 0.5d0*rhor*ur**2
c
      return
      end
