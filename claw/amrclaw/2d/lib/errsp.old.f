c
c ----------------------------------------------------------
c
      subroutine errsp(rect,sperr,mitot,mjtot,nvar,mptr,ng,
     &                 eprint, outunit)

c
c ::::::::::::::::::::: ERRSP ::::::::::::::::::::::::::::::::::
c estimate spatial only component of the error
c rect   = grid values including ghost cells
c sperr  = user computed fn. (density gradient here)
c           (if sperr(i,j) > tolsp cell will be flagged)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit double precision (a-h, o-z)

c      include  "call.i"

      dimension   rect(mitot,mjtot,nvar)
      dimension  sperr(mitot,mjtot)
      integer    outunit
      logical    eprint

      drhomax = 0.d0
      imax    = 0
      jmax    = 0

      do 20 j = ng+1, mjtot-ng
      do 10 i = ng+1, mitot-ng

c        # look for energy jumps too!
         drhoi = dabs(rect(i+1,j,1) - rect(i-1,j,1))
         drhoj = dabs(rect(i,j+1,1) - rect(i,j-1,1))
         drho  = dmax1(drhoi, drhoj)
         sperr(i,j) = drho
         if (drho .gt. drhomax) then
           drhomax = drho
           imax = i
           jmax = j
         endif

 10   continue
 20   continue

      if (eprint) then
          write(outunit,*) ' maximum drho for grid ',mptr, 
     .                     ' was ', drhomax
          write(outunit,*)    ' at location ',imax,jmax
      endif

      return
      end
