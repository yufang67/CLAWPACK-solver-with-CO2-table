c
c --------------------------------------------------------------
c
      subroutine prefil4(level,nvar,
     1                   valbig,aux,naux,time,mitot,mjtot,
     1                   nrowst,ncolst,
     3                   ilo,ihi,jlo,jhi)

c
c  6/21/05: added another level of filpatch and prefil since a 4th level
c  may be needed with 2 ghost cells and refinement by 2 and centered 
c  interpolation for ghost values.  We put in files with filpatch3 (prefil3) 
c  so the Makefiles do not change.

      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c
c  :::::::::::::: PREFILP :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to fill the
c     piece of the boundary. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     filpatch to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     xl, xr, yb, yt = the location in physical space of
c     corners of a patch.
c     ilo,ihi,jlo,jhi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c
c     # will divide patch into 9 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
        ist(1) = ilo
        ist(2) = 0
        ist(3) = iregsz(level)
        iend(1) = -1
        iend(2) = iregsz(level)-1
        iend(3) = ihi
        jst(1) = jlo
        jst(2) = 0
        jst(3) = jregsz(level)
        jend(1) = -1
        jend(2) = jregsz(level)-1
        jend(3) = jhi
        ishift(1) = iregsz(level)
        ishift(2) = 0
        ishift(3) = -iregsz(level)
        jshift(1) = jregsz(level)
        jshift(2) = 0
        jshift(3) = -jregsz(level)

       do 20 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
       do 10 j = 1, 3
          j1 = max(jlo,  jst(j))
          j2 = min(jhi, jend(j))


          if ((i1 .le. i2) .and. (j1 .le. j2)) then
            iputst = (i1 - ilo) + nrowst
            jputst = (j1 - jlo) + ncolst
            call filpatch4(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                     iputst,jputst,
     2                    i1+ishift(i),i2+ishift(i),
     3                    j1+jshift(j),j2+jshift(j))
          endif

 10    continue
 20    continue
      
     


      return
      end
c
c --------------------------------------------------------------
c
      subroutine prefil3(level,nvar,
     1                   valbig,aux,naux,time,mitot,mjtot,
     1                   nrowst,ncolst,
     3                   ilo,ihi,jlo,jhi)

c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c
c  :::::::::::::: PREFILP :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to fill the
c     piece of the boundary. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     filpatch to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     xl, xr, yb, yt = the location in physical space of
c     corners of a patch.
c     ilo,ihi,jlo,jhi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c
c     # will divide patch into 9 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
        ist(1) = ilo
        ist(2) = 0
        ist(3) = iregsz(level)
        iend(1) = -1
        iend(2) = iregsz(level)-1
        iend(3) = ihi
        jst(1) = jlo
        jst(2) = 0
        jst(3) = jregsz(level)
        jend(1) = -1
        jend(2) = jregsz(level)-1
        jend(3) = jhi
        ishift(1) = iregsz(level)
        ishift(2) = 0
        ishift(3) = -iregsz(level)
        jshift(1) = jregsz(level)
        jshift(2) = 0
        jshift(3) = -jregsz(level)

       do 20 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
       do 10 j = 1, 3
          j1 = max(jlo,  jst(j))
          j2 = min(jhi, jend(j))


          if ((i1 .le. i2) .and. (j1 .le. j2)) then
            iputst = (i1 - ilo) + nrowst
            jputst = (j1 - jlo) + ncolst
            call filpatch3(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                     iputst,jputst,
     2                    i1+ishift(i),i2+ishift(i),
     3                    j1+jshift(j),j2+jshift(j))
          endif

 10    continue
 20    continue
      
     


      return
      end
