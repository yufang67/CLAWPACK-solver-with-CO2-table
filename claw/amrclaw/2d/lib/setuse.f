c
c ------------------------------------------------------------------
c
       subroutine setuse(listbc,maxsp,ispot,mkid,
     1                   ilo,ihi,jlo,jhi,
     2                   iclo,ichi,jclo,jchi,nghost)
c
c ::::::::::::::::::::::::: SETUSE ::::::::::::::::::::::::::::::::
c
c set up boundary list for coarse grid, to be used by fluxsv. 
c loop around boundary of c fine grids to do this.  each entry has
c     i, j, side #, fine grid #, loc in fine grid list for fluxes.
c  for example, side 1 of fine grid fixes side 3 of coarse grid,
c  so coarse grid list will store the # 3.
c  wrt coarse grid, the sides are:
c              2
c           1     3            that is, right edge of a coarse cell = 3
c              4                         top  edge of a coarse cell = 2
c
c  # lkid is the index into the fine grid's saved fluxes.
c  # the fine grid will save all its fluxes all around its
c  # perimeter. lkid tells where the coarse grid should
c  # taking them from.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      implicit double precision (a-h,o-z)
      dimension listbc(5,maxsp)

      ibc = ispot
      ist  = iclo - 1
      iend = ichi + 1
      jst  = jclo - 1
      jend = jchi + 1
c
c  left side
c
       if (ist .lt. ilo) go to 20
       lkid     = max(jlo,jclo) - jclo + 1
       do 10 j  = max(jlo,jclo), min(jhi,jchi)
          ispot              = ispot + 1
          listbc(1,ispot)    = ist-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = 3
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid               = lkid + 1
 10    continue
c
c   top side
c
 20    if (jend .gt. jhi) go to 40
       lkid       = (jchi-jclo+1) + max(ilo,iclo)-iclo + 1
       do 30 i    = max(ilo,iclo), min(ihi,ichi)
          ispot              = ispot + 1
          listbc(1,ispot)    = i-ilo+nghost+1
          listbc(2,ispot)    = jend-jlo+nghost+1
          listbc(3,ispot)    = 4
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid               = lkid + 1
 30    continue
c
c  right side (numbered from bottom to top, so not continuous)
c
 40    if (iend .gt. ihi) go to 60
       lkid     = (ichi-iclo+1)+(jchi-jclo+1)
     .               + max(jlo,jclo) - jclo + 1
       do 50 j  = max(jlo,jclo), min(jhi,jchi)
          ispot              = ispot + 1
          listbc(1,ispot)    = iend-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = 1
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid   = lkid + 1
 50    continue
c
c  bottom side (numbered left to right, so not continuous)
c
 60    if (jst .lt. jlo) go to 80
       lkid   =  2*(jchi-jclo+1) + (ichi-iclo+1) + max(ilo,iclo)-iclo+1
       do 70 i          = max(ilo,iclo), min(ihi,ichi)
          ispot              = ispot + 1
          listbc(1,ispot)    = i-ilo+nghost+1
          listbc(2,ispot)    = jst-jlo+nghost+1
          listbc(3,ispot)    = 2
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid   = lkid + 1
 70    continue
c
 80    continue
       return
       end
