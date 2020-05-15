      program driver
c
c  Generic driver routine for claw2
c
c  Author: Randall J. LeVeque
c  Version of March, 1999 --  CLAWPACK Version 4.0
c
c
      implicit double precision (a-h,o-z)

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c
      parameter (maxmx =   1000)
      parameter (maxmy =   1000)
      parameter (mwork =  9100599)

      parameter (mbc = 2)
      parameter (meqn = 4)
      parameter (mwaves = 3)
      parameter (maux = 7)
c       # NOTE: if maux>0 you must declare aux properly below!

      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

c      dimension  aux(1)   !# dummy variable since no aux arrays used
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
c     dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      dimension mthlim(mwaves)
      dimension work(mwork)
c
      call claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,aux)

      stop 
      end
