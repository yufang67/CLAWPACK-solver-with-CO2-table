c
c ------------------------------------------------
c
      subroutine drawrg(time,lpar,mkid1,nclust,numptc,nxypts,
     *                  badpts)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension badpts(2,nxypts), numptc(nclust)

c :::::::::::::::::::::: DRAWRG ::::::::::::::::::::::::::::::;
c
c drawrg: output existing grids at level lpar
c                flagged points at level lpar
c                clusters
c                new and old subgrids - start with grid mkid1 (not a
c                     levelptr.)
c  make sure nclust <> 0 at this point. mkid1 always exists if
c  there is at least 1 flagged pt., but there might not be any
c  old fine grids.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      if (nclust .le. 0) return
c
      call basic (time, lpar, lpar+1)
c
      write(pltunit1,100) nclust, lpar, nxypts
100   format(10h*CLUSTERS ,3i10)
      write(pltunit1,105) iregsz, jregsz
105   format(6i10)
      write(pltunit1,101) (numptc(i),i=1, nclust)
101   format(10i6)
      write(pltunit1,102) ((badpts(i,j),i=1,2),j=1,nxypts)
102   format(4e12.4)
c
c output subgrid info. - since nclust > 0, we know there is at least
c  one. then we can chain through levelptr to get rest.
c
      mkid = mkid1
 10   if (mkid .eq. 0) go to 20
          write(pltunit1,103) mkid, (node(i,mkid),i=1,nsize)
103       format(10i6)
          write(pltunit1,104) (rnode(i,mkid), i = 1,rsize)
104       format(5e12.4)
          mkid = node(levelptr, mkid)
      go to 10
 20   continue
c
      return
      end
