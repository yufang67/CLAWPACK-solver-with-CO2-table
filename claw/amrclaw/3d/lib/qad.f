c
c -------------------------------------------------------------
c
       subroutine qad(valbig,mitot,mjtot,mktot,nvar,
     .                svdflx,qc1d,lenbc,lratiox,lratioy,lratioz,
     .                hx,hy,hz,maux,aux,auxc1d,delt,mptr)

       implicit double precision (a-h, o-z)

       include "call.i"

       logical qprint

       dimension valbig(mitot,mjtot,mktot,nvar)
       dimension qc1d(lenbc,nvar)
       dimension svdflx(nvar,lenbc)
       dimension aux(mitot,mjtot,mktot,maux)
       dimension auxc1d(lenbc,maux)

c
c ::::::::::::::::::::::::::: QAD ::::::::::::::::::::::::::::::::::
c  solve RP between ghost cell value on fine grid and coarse grid
c  value that ghost cell overlaps. The resulting fluctuations
c  are added in to coarse grid value, as a conservation fixup.
c  Done each fine grid time step. If source terms are present, the
c  coarse grid value is advanced by source terms each fine time step too.

c Side 1 is the left   side of the fine grid patch.   i = min
c Side 2 is the rear   side of the fine grid patch.   j = max
c Side 3 is the right  side of the fine grid patch.   i = max
c Side 4 is the front  side of the fine grid patch.   j = min
c Side 5 is the bottom side of the fine grid patch.   k = min
c Side 6 is the top    side of the fine grid patch.   k = max
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c      # local storage
c      # note that dimension here are bigger than dimensions used
c      # in rp2, but shouldn't matter since wave is not used in qad
c      # and for other arrays it is only the last parameter that is wrong
c      #  ok as long as meqn, mwaves < maxvar

       parameter (max1dp1 = max1d+1)
c      parameter (max1dp1 = (max1d+1)**2+1)
       dimension ql(max1dp1,maxvar),    qr(max1dp1,maxvar)
       dimension wave(max1dp1,maxvar,maxvar),  s(max1dp1,maxvar)
       dimension amdq(max1dp1,maxvar),  apdq(max1dp1,maxvar)
       dimension auxl(max1dp1,maxaux),  auxr(max1dp1,maxaux)

       data qprint/.false./
c
c      aux is auxiliary array with user parameters needed in Riemann solvers
c          on fine grid corresponding to valbig
c      auxc1d is coarse grid stuff from around boundary, same format as qc1d
c      auxl, auxr are work arrays needed to pass stuff to rpn2
c      maux is the number of aux variables, which may be zero.
c
c      nr, nc, nf are the number of rows, columns, files
c


       if (qprint) write(dbugunit,*)" working on grid ",mptr
       tgrid = rnode(timemult, mptr)
       nf = mktot-2*nghost
       nc = mjtot-2*nghost
       nr = mitot-2*nghost

       ncrsek = (mktot-2*nghost)/lratioz
       ncrsej = (mjtot-2*nghost)/lratioy
       ncrsei = (mitot-2*nghost)/lratiox
c
c--------
c  side 1
c--------
c
      index0  = 0
      influx0 = 0
      do 199 k = nghost+1, mktot-nghost
      kc     = 1 + (k  -(nghost+1))/lratioz
      kcnext = 1 + (k+1-(nghost+1))/lratioz
      indaux = 1
      do 10 j = nghost+1, mjtot-nghost
      indaux = indaux + 1
      if (maux.gt.0) then
         do 5 ma = 1,maux
	    if (auxtype(ma).eq."xleft") then
c               # Assuming velocity at left-face, this fix
c               # preserves conservation in incompressible flow:
   	        auxl(indaux,ma) = aux(nghost+1,j,k,ma)
	    else
c               # Normal case -- we set the aux arrays
c               # from the cell corresponding  to q
	        auxl(indaux,ma) = aux(nghost,j,k,ma)
	    endif
    5    continue
      endif
      do 10 ivar = 1, nvar
         ql(indaux,ivar) = valbig(nghost,j,k,ivar)
   10 continue

      lind  = 0
      index = index0
      do 22 jc = 1, ncrsej
        index = index + 1
        do 23 l = 1, lratioy
          lind = lind + 1
          if (maux.gt.0) then
            do 24 ma=1,maux
              auxr(lind,ma) = auxc1d(index,ma)
   24       continue
          endif
          do 25 ivar = 1, nvar
            qr(lind,ivar) = qc1d(index,ivar)
   25     continue
   23   continue
   22 continue
      if (kcnext .ne. kc) index0 = index

       if (qprint) then
	 write(dbugunit,*) 'side 1, ql and qr:'
	 do i=2,nc
	    write(dbugunit,4101) i,qr(i-1,1),ql(i,1),auxr(i-1,1),auxl(i,1)
 	  enddo
 4101      format(i3,4e16.6)
       endif

       call rpn3(1,max1dp1-2*nghost,nvar,mwaves,nghost,nc+1-2*nghost,
     .              ql,qr,auxl,auxr,maux,wave,s,amdq,apdq)
c
c we have the wave. for side 1 add into sdflxm
c
      kfine   = (kc-1)*lratioz
      influx  = influx0
      do 50 j = 1, nc/lratioy
         influx  = influx + 1
         jfine   = (j-1)*lratioy
*        indexa  = (kfine+m-1)*nc + jfine + 1
         indexa  =                  jfine + 1
         do 60 ivar = 1, nvar
            do 70 l = 1, lratioy
               svdflx(ivar,influx) = svdflx(ivar,influx)
     .                       + amdq(indexa+l,ivar) * hy * hz * delt
     .                       + apdq(indexa+l,ivar) * hy * hz * delt
   70       continue
   60    continue
   50 continue
      if (kcnext .ne. kc) influx0 = influx

  199 continue

c--------
c  side 2
c--------
c
      do 299 k = nghost+1, mktot-nghost
      kc     = 1 + (k  -(nghost+1))/lratioz
      kcnext = 1 + (k+1-(nghost+1))/lratioz
      indaux = 0
      do 210 i = nghost+1, mitot-nghost
      indaux = indaux + 1
      if (maux.gt.0) then
         do 205 ma = 1,maux
            auxr(indaux,ma) = aux(i,mjtot-nghost+1,k,ma)
  205    continue
      endif
      do 210 ivar = 1, nvar
         qr(indaux,ivar) = valbig(i,mjtot-nghost+1,k,ivar)
  210 continue

      lind  = 0
      index = index0
      do 222 ic = 1, ncrsei
        index = index + 1
        do 223 l = 1, lratiox
          lind = lind + 1
          if (maux.gt.0) then
            do 224 ma=1,maux
              if (auxtype(ma).eq."yleft") then
c                 # Assuming velocity at rear-face, this fix
c                 # preserves conservation in incompressible flow:
                  kfine           = (kc-1)*lratioz + nghost + l
                  ifine           = (ic-1)*lratiox + nghost + l
                  auxl(lind+1,ma) = aux(ifine,mjtot-nghost+1,kfine,ma)
              else
                  auxl(lind+1,ma) = auxc1d(index,ma)
              endif
  224       continue
          endif
          do 225 ivar = 1, nvar
            ql(lind+1,ivar) = qc1d(index,ivar)
  225     continue
  223   continue
  222 continue
      if (kcnext .ne. kc) index0 = index

       if (qprint) then
	 write(dbugunit,*) 'side 2, ql and qr:'
	 do i=2,nr
	    write(dbugunit,4101) i,qr(i-1,1),ql(i,1),auxr(i-1,1),auxl(i,1)
	    enddo
       endif

       call rpn3(2,max1dp1-2*nghost,nvar,mwaves,nghost,nr+1-2*nghost,
     .              ql,qr,auxl,auxr,maux,wave,s,amdq,apdq)
c
c we have the wave. for side 2. add into sdflxp
c
      kfine   = (kc-1)*lratioz
      influx  = influx0
      do 250 i = 1, nr/lratiox
         influx  = influx + 1
         ifine   = (i-1)*lratiox
*        indexa  = (kfine+m-1)*nc + ifine + 1
         indexa  =                  ifine + 1
         do 260 ivar = 1, nvar
            do 270 l = 1, lratiox
               svdflx(ivar,influx) = svdflx(ivar,influx)
     .                       - amdq(indexa+l,ivar) * hx * hz * delt
     .                       - apdq(indexa+l,ivar) * hx * hz * delt
  270       continue
  260    continue
  250 continue
      if (kcnext .ne. kc) influx0 = influx

  299 continue

c--------
c  side 3
c--------
c
      do 399 k = nghost+1, mktot-nghost
      kc     = 1 + (k  -(nghost+1))/lratioz
      kcnext = 1 + (k+1-(nghost+1))/lratioz
      indaux = 0
      do 310 j = nghost+1, mjtot-nghost
      indaux = indaux + 1
      if (maux.gt.0) then
         do 305 ma = 1,maux
            auxr(indaux,ma) = aux(mitot-nghost+1,j,k,ma)
  305    continue
      endif
      do 310 ivar = 1, nvar
         qr(indaux,ivar) = valbig(mitot-nghost+1,j,k,ivar)
  310 continue

      lind  = 0
      index = index0
      do 322 jc = 1, ncrsej
        index = index + 1
        do 323 l = 1, lratioy
          lind = lind + 1
          if (maux.gt.0) then
            do 324 ma=1,maux
              if (auxtype(ma).eq."xleft") then
c               # Assuming velocity at left-face, this fix
c               # preserves conservation in incompressible flow:
                kfine           = (kc-1)*lratioz + nghost + l
                jfine           = (jc-1)*lratioy + nghost + l
                auxl(lind+1,ma) = aux(mitot-nghost+1,jfine,kfine,ma)
              else
                auxl(lind+1,ma) = auxc1d(index,ma)
              endif
  324       continue
          endif
          do 325 ivar = 1, nvar
             ql(lind+1,ivar) = qc1d(index,ivar)
  325     continue
  323   continue
  322 continue
      if (kcnext .ne. kc) index0 = index

       if (qprint) then
	 write(dbugunit,*) 'side 3, ql and qr:'
	 do i=2,nc
	    write(dbugunit,4101) i,qr(i-1,1),ql(i,1),auxr(i-1,1),auxl(i,1)
	    enddo
       endif

       call rpn3(1,max1dp1-2*nghost,nvar,mwaves,nghost,nc+1-2*nghost,
     .              ql,qr,auxl,auxr,maux,wave,s,amdq,apdq)
c
c we have the wave. for side 3 add into sdflxp
C
      kfine   = (kc-1)*lratioz
      influx = influx0
      do 350 j = 1, nc/lratioy
         influx  = influx + 1
         jfine   = (j-1)*lratioy
*        indexa  = (kfine+m-1)*nc + jfine + 1
         indexa  =                  jfine + 1
         do 360 ivar = 1, nvar
            do 370 l = 1, lratioy
               svdflx(ivar,influx) = svdflx(ivar,influx)
     .                       - amdq(indexa+l,ivar) * hy * hz * delt
     .                       - apdq(indexa+l,ivar) * hy * hz * delt
  370       continue
  360    continue
  350 continue
      if (kcnext .ne. kc) influx0 = influx

  399 continue

c--------
c  side 4
c--------
c
      do 499 k = nghost+1, mktot-nghost
      kc     = 1 + (k  -(nghost+1))/lratioz
      kcnext = 1 + (k+1-(nghost+1))/lratioz
      indaux = 1
      do 410 i = nghost+1, mitot-nghost
      indaux = indaux + 1
      if (maux.gt.0) then
         do 405 ma = 1,maux
            if (auxtype(ma).eq."yleft") then
c               # Assuming velocity at rear-face, this fix
c               # preserves conservation in incompressible flow:
                auxl(indaux,ma) = aux(i,nghost+1,k,ma)
            else
                auxl(indaux,ma) = aux(i,nghost,k,ma)
            endif
  405    continue
      endif
      do 410 ivar = 1, nvar
         ql(indaux,ivar) = valbig(i,nghost,k,ivar)
  410 continue

      lind  = 0
      index = index0
      do 422 ic = 1, ncrsei
        index = index + 1
        do 423 l = 1, lratiox
          lind = lind + 1
          if (maux.gt.0) then
            do 424 ma=1,maux
              auxr(lind,ma) = auxc1d(index,ma)
  424       continue
          endif
          do 425 ivar = 1, nvar
            qr(lind,ivar) = qc1d(index,ivar)
  425     continue
  423   continue
  422 continue
      if (kcnext .ne. kc) index0 = index

       if (qprint) then
	 write(dbugunit,*) 'side 4, ql and qr:'
	 do i=2,nr
	    write(dbugunit,4101) i,qr(i-1,1),ql(i,1),auxr(i-1,1),auxl(i,1)
	    enddo
       endif

       call rpn3(2,max1dp1-2*nghost,nvar,mwaves,nghost,nr+1-2*nghost,
     .              ql,qr,auxl,auxr,maux,wave,s,amdq,apdq)
c
c we have the wave. for side 4. add into sdflxm
c
      kfine   = (kc-1)*lratioz
      influx  = influx0
      do 450 i = 1, nr/lratiox
         influx  = influx + 1
         ifine   = (i-1)*lratiox
*        indexa  = (kfine+m-1)*nc + ifine + 1
         indexa  =                  ifine + 1
         do 460 ivar = 1, nvar
            do 470 l = 1, lratiox
               svdflx(ivar,influx) = svdflx(ivar,influx)
     .                       + amdq(indexa+l,ivar) * hx * hz * delt
     .                       + apdq(indexa+l,ivar) * hx * hz * delt
  470       continue
  460    continue
  450 continue
      if (kcnext .ne. kc) influx0 = influx

  499 continue

c--------
c  side 5
c--------
c
      do 599 j = nghost+1, mjtot-nghost
      jc     = 1 + (j  -(nghost+1))/lratioy
      jcnext = 1 + (j+1-(nghost+1))/lratioy
      indaux = 1
      do 510 i = nghost+1, mitot-nghost
      indaux = indaux + 1
      if (maux.gt.0) then
         do 505 ma = 1,maux
            if (auxtype(ma).eq."zleft") then
c               # Assuming velocity at bottom-face, this fix
c               # preserves conservation in incompressible flow:
                auxl(indaux,ma) = aux(i,j,nghost+1,ma)
            else
c               # Normal case -- we set the aux arrays
c               # from the cell corresponding  to q
                auxl(indaux,ma) = aux(i,j,nghost,ma)
            endif
  505       continue
      endif
      do 510 ivar = 1, nvar
         ql(indaux,ivar) = valbig(i,j,nghost,ivar)
  510 continue

      lind  = 0
      index = index0
      do 522 ic = 1, ncrsei
        index = index + 1
        do 523 l = 1, lratiox
          lind = lind + 1
          if (maux.gt.0) then
            do 524 ma=1,maux
              auxr(lind,ma) = auxc1d(index,ma)
  524       continue
          endif
          do 525 ivar = 1, nvar
            qr(lind,ivar) = qc1d(index,ivar)
  525     continue
  523   continue
  522 continue
      if (jcnext .ne. jc) index0 = index

       if (qprint) then
         write(dbugunit,*) 'side 5, ql and qr:'
         do i=2,nc
	    write(dbugunit,4101) i,qr(i-1,1),ql(i,1),auxr(i-1,1),auxl(i,1)
          enddo
       endif

       call rpn3(3,max1dp1-2*nghost,nvar,mwaves,nghost,nr+1-2*nghost,
     .              ql,qr,auxl,auxr,maux,wave,s,amdq,apdq)
c
c we have the wave. for side 5 add into sdflxm
c
      jfine   = (jc-1)*lratioy
      influx  = influx0
      do 550 i = 1, nr/lratiox
         influx  = influx + 1
         ifine   = (i-1)*lratiox
*        indexa  = (jfine+m-1)*nr + ifine + 1
         indexa  =                  ifine + 1
         do 560 ivar = 1, nvar
            do 570 l = 1, lratiox
               svdflx(ivar,influx) = svdflx(ivar,influx)
     .                       + amdq(indexa+l,ivar) * hx * hy * delt
     .                       + apdq(indexa+l,ivar) * hx * hy * delt
  570       continue
  560    continue
  550 continue
      if (jcnext .ne. jc) influx0 = influx

  599 continue

c--------
c  side 6
c--------
c
      do 699 j = nghost+1, mjtot-nghost
      jc     = 1 + (j  -(nghost+1))/lratioy
      jcnext = 1 + (j+1-(nghost+1))/lratioy
      indaux = 0
      do 610 i = nghost+1, mitot-nghost
      indaux = indaux + 1
      if (maux.gt.0) then
         do 605 ma = 1,maux
            auxr(indaux,ma) = aux(i,j,mktot-nghost+1,ma)
  605    continue
      endif
      do 610 ivar = 1, nvar
         qr(indaux,ivar) = valbig(i,j,mktot-nghost+1,ivar)
  610 continue

      lind   = 0
      index = index0
      do 622 ic = 1, ncrsei
        index = index + 1
        do 623 l = 1, lratiox
          lind = lind + 1
          if (maux.gt.0) then
            do 624 ma=1,maux
              if (auxtype(ma).eq."zleft") then
c                 # Assuming velocity at bottom-face, this fix
c                 # preserves conservation in incompressible flow:
                  jfine = (jc-1)*lratioy + nghost + l
                  ifine = (ic-1)*lratiox + nghost + l
                  auxl(lind+1,ma) = aux(ifine,jfine,mktot-nghost+1,ma)
              else
                  auxl(lind+1,ma) = auxc1d(index,ma)
              endif
  624       continue
          endif
          do 625 ivar = 1, nvar
            ql(lind+1,ivar) = qc1d(index,ivar)
  625     continue
  623   continue
  622 continue
      if (jcnext .ne. jc) index0 = index

       if (qprint) then
         write(dbugunit,*) 'side 6, ql and qr:'
         do i=2,nr
	    write(dbugunit,4101) i,qr(i-1,1),ql(i,1),auxr(i-1,1),auxl(i,1)
            enddo
       endif

       call rpn3(3,max1dp1-2*nghost,nvar,mwaves,nghost,nr+1-2*nghost,
     .              ql,qr,auxl,auxr,maux,wave,s,amdq,apdq)
c
c we have the wave. for side 6 add into sdflxp
C
      jfine   = (jc-1)*lratioy
      influx = influx0
      do 650 i = 1, nr/lratiox
         influx  = influx + 1
         ifine   = (i-1)*lratiox
*        indexa  = (jfine+m-1)*nr + ifine + 1
         indexa  =                  ifine + 1
         do 660 ivar = 1, nvar
            do 670 l = 1, lratiox
               svdflx(ivar,influx) = svdflx(ivar,influx)
     .                       - amdq(indexa+l,ivar) * hx * hy * delt
     .                       - apdq(indexa+l,ivar) * hx * hy * delt
  670       continue
  660    continue
  650 continue
      if (jcnext .ne. jc) influx0 = influx

  699 continue

c      # for source terms:
       if (method(5) .ne. 0) then
c	   call src3(qc1d,lenbc,1,nvar,auxc1d,maux,tgrid,delt)
           call src1d(nvar,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
       endif

       return
       end
