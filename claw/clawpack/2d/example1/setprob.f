      subroutine setprob
      implicit double precision (a-h,o-z)
      common /comrp/ ubar,vbar
c
c     # Set the velocity for scalar advection
c     # These values are passed to the Riemann solvers rpn2.f and rpt2.f
c     # in a common block
c

      open(unit=7,file='setprob.data',status='old',form='formatted')

      read(7,*) ubar
      read(7,*) vbar

      return
      end
