c
c ----------------------------------------------------------
c
      subroutine errsp(q,sperr,mitot,mjtot,nvar,nghost,
     &                 tolsp,time,GOODFLAG,BADFLAG)

c **********************
c ****** OBSOLETE ******  No longer used in amrclaw.
c **********************
c
c ::::::::::::::::::::: ERRSP ::::::::::::::::::::::::::::::::::
c user routine to control flagging of points based on
c gradients or other simple criteria
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit double precision (a-h, o-z)

      dimension   q(mitot,mjtot,nvar)
      dimension   sperr(mitot,mjtot)
c
c     # This routine is obsolete and no longer used.
c
      write(6,*) '*** The new version of AMRCLAW requires some changes'
      write(6,*) 'to Makefile:'
      write(6,*) '   Remove errsp.f (no longer used)'
      write(6,*) '   Add the following library routines to Makefile:'
      write(6,*) '      $(CLAW)/amrclaw/2d/lib/bufnst.f'
      write(6,*) '      $(CLAW)/amrclaw/2d/lib/spest.f'
      write(6,*) '      $(CLAW)/amrclaw/2d/lib/flag2refine.f'
      write(6,*) '      $(CLAW)/amrclaw/2d/lib/allowflag.f'
      write(6,*) 'See the documentation in flag2refine.f and '
      write(6,*) '    allowflag.f for more information.'
      stop

      end
