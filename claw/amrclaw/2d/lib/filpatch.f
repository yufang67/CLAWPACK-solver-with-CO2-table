c
c ---------------------------------------------------------------
c
        subroutine filpatch(level,nvar,valbig,aux,naux,
     1                      time,mitot,mjtot,
     2                      nrowst,ncolst,ilo,ihi,jlo,jhi)

c :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
c
c  fill the portion of valbig from rows  nrowst
c                             and  cols  ncolst
c  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
c  vals are needed at time time, and level level,
c
c  first fill with  values obtainable from the level level
c  grids. if any left unfilled, then enlarge remaining rectangle of
c  unfilled values by 1 (for later linear interp), and recusively
c  obtain the remaining values from  coarser levels.
c
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;

      implicit double precision (a-h,o-z)

      include  "call.i"

      logical   set, sticksout
      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)

      iadflag(i,j)    =  locuse + i-1+(j-1)*nrowp
      ivalc(i,j,ivar) =  loccrse + (i - 1) + nrowc*(j - 1)
     &                     + nrowc*ncolc*(ivar-1)
      sticksout(iplo,iphi,jplo,jphi)  =
     &            (iplo .lt. 0 .or. jplo .lt. 0 .or.
     &             iphi .ge. iregsz(levc) .or. jphi .ge. jregsz(levc))
c
c We begin by filling values for grids at level level. If all values can be
c filled in this way, we return;

        hxf     = hxposs(level)
        hyf     = hyposs(level)
        xlp     = xlower + ilo*hxf
        xrp     = xlower + (ihi+1)*hxf
        ybp     = ylower + jlo*hyf
        ytp     = ylower + (jhi+1)*hyf
        nrowp   = ihi - ilo + 1
        ncolp   = jhi - jlo + 1
        locuse  = igetsp(nrowp*ncolp)

c 2/27/02 : added naux to intfil call.
        call intfil
     &  (valbig,mitot,mjtot,time,locuse,nrowst,ncolst,
     &   ilo,ihi,jlo,jhi,level,nvar,naux)
c
c Trimbd returns set = true if all of the entries are filled (=1.).
c set = false, otherwise. If set = true, then no other levels are
c are required to interpolate, and we return.
c
c Note that the used array is filled entirely in intfil, i.e. the
c marking done there also takes into account the points filled by
c the boundary conditions. bc2amr will be called later (from bound), after
c all 4 boundary pieces filled.

        call trimbd(alloc(locuse),nrowp,ncolp,set,il,ir,jb,jt)

        if (set) then
           call reclam(locuse,nrowp*ncolp)
           return
        else if (level .eq. 1) then
           write(outunit,*)" error in filpatch - level 1 not set"
           write(outunit,900) nrowst,ncolst
           write(*,*)" error in filpatch - level 1 not set"
           write(*,900) nrowst,ncolst
900        format("starting at row: ",i4," col ",i4)
           stop
        endif

c set = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c recursive formulation for interpolating.
c IS THIS TRUE ? - the fine grid patch remaining unfilled is always
c anchored to a coarse cell.

        levc = level - 1
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)

        isl  = il + ilo - 1
        isr  = ir + ilo - 1
        jsb  = jb + jlo - 1
        jst  = jt + jlo - 1
c
c       coarsen
        lratiox = intratx(levc)
        lratioy = intraty(levc)

        iplo   = (isl-lratiox  +nghost*lratiox)/lratiox - nghost
        jplo   = (jsb-lratioy  +nghost*lratioy)/lratioy - nghost
        iphi   = (isr+lratiox  )/lratiox
        jphi   = (jst+lratioy  )/lratioy


        xlc  =  xlower + iplo*hxc
        ybc  =  ylower + jplo*hyc
        xrc  =  xlower + (iphi+1)*hxc
        ytc  =  ylower + (jphi+1)*hyc

        nrowc   =  iphi - iplo + 1
        ncolc   =  jphi - jplo + 1
        ntot    = nrowc*ncolc*(nvar+naux)
        loccrse = igetsp(ntot)
        locauxc = loccrse + nrowc*ncolc*nvar
        if (naux.gt.0) then
              maxmx = nrowc - 2*nghost
              mx = maxmx
              maxmy = ncolc - 2*nghost
              my = maxmy
              xl = xlc + nghost*hxc
              yb = ybc + nghost*hyc
              call setaux(maxmx,maxmy,nghost,mx,my,xl,yb,hxc,hyc,
     &                    naux,alloc(locauxc))
        endif


        if ((xperdom .or. yperdom) .and.
     &       sticksout(iplo,iphi,jplo,jphi)) then
             call prefil2(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
     1                    time,nrowc,ncolc,1,1,
     2                    iplo,iphi,jplo,jphi)
        else
          call filpatch2(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
     1                   time,nrowc,ncolc,1,1,
     2                   iplo,iphi,jplo,jphi)
        endif


c       interpolate back up

20      continue

        do 100 iff = 1,nrowp
          ic = 2 + (iff - (isl - ilo) - 1)/lratiox
          eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)

        do 100 jf  = 1,ncolp
          jc = 2 + (jf -(jsb-jlo)-1)/lratioy
          eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)

           flag = alloc(iadflag(iff,jf))
           if (flag .eq. 0.0) then

c               xif = xlp + (.5 + float(iff-1))*hxf
c               yjf = ybp + (.5 + float(jf -1))*hyf

c               ic=idint((xif-xlc+.5*hxc)/hxc)
c               jc=idint((yjf-ybc+.5*hyc)/hyc)

c               xc = xlc + (float(ic) - .5)*hxc
c               yc = ybc + (float(jc) - .5)*hyc

c               eta1 = (xif - xc)/hxc
c               eta2 = (yjf - yc)/hyc

                do 101 ivar = 1,nvar

c                  valc00 = alloc(ivalc(ic,jc,ivar))
c                  valc10 = alloc(ivalc(ic+1,jc,ivar))
c                  valc01 = alloc(ivalc(ic,jc+1,ivar))
c                  valc11 = alloc(ivalc(ic+1,jc+1,ivar))

                   valp10 = alloc(ivalc(ic+1,jc,ivar))
                   valm10 = alloc(ivalc(ic-1,jc,ivar))
                   valc   = alloc(ivalc(ic  ,jc,ivar))
                   valp01 = alloc(ivalc(ic  ,jc+1,ivar))
                   valm01 = alloc(ivalc(ic  ,jc-1,ivar))

                   dupc = valp10 - valc
                   dumc = valc   - valm10
                   ducc = valp10 - valm10
                   du   = dmin1(dabs(dupc),dabs(dumc))
                   du   = dmin1(2.d0*du,.5d0*dabs(ducc))
                   fu = dmax1(0.d0,dsign(1.d0,dupc*dumc))

                   dvpc = valp01 - valc
                   dvmc = valc   - valm01
                   dvcc = valp01 - valm01
                   dv   = dmin1(dabs(dvpc),dabs(dvmc))
                   dv   = dmin1(2.d0*dv,.5d0*dabs(dvcc))
                   fv = dmax1(0.d0,dsign(1.d0,dvpc*dvmc))

                   valint = valc + eta1*du*dsign(1.d0,ducc)*fu
     .                           + eta2*dv*dsign(1.d0,dvcc)*fv

c                  valint = (1. - eta2)*
c    &               ((1. - eta1)*valc00 + eta1*valc10)
c    &               + eta2*((1. - eta1)*valc01 + eta1*valc11)

                   valbig(iff+nrowst-1,jf+ncolst-1,ivar) = valint

101             continue

           endif

100     continue

        call reclam(loccrse,ntot)

        call reclam(locuse,nrowp*ncolp)

        return
        end
