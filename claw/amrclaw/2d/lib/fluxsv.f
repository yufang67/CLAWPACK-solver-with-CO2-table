c
c ----------------------------------------------------------
c
      subroutine fluxsv(mptr,xfluxm,xfluxp,yfluxm,yfluxp,listbc,
     1                  ndimx,ndimy,nvar,maxsp,dtc,hx,hy)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension xfluxp(ndimx,ndimy,nvar), yfluxp(ndimx,ndimy,nvar)
      dimension xfluxm(ndimx,ndimy,nvar), yfluxm(ndimx,ndimy,nvar)
      dimension listbc(5,maxsp)
c
c :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
c
c  coarse grids should save their fluxes in cells adjacent to
c  their nested fine grids, for later conservation fixing.
c  listbc holds info for where to save which fluxes.
c  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
c
c :::::::::::::::::::::::::::::;:::::::::::::::::::::::
 
 
      ispot   = 1
      level   = node(nestlevel,mptr)

 10      if (listbc(1,ispot).eq.0) go to 99          
c
         mkid     = listbc(4,ispot)
         intopl   = listbc(5,ispot)
         nx       = node(ndihi,mkid) - node(ndilo,mkid) + 1
         ny       = node(ndjhi,mkid) - node(ndjlo,mkid) + 1
         kidlst   = node(ffluxptr,mkid)
         i        = listbc(1,ispot)
         j        = listbc(2,ispot)
         inlist   = kidlst + nvar*(intopl-1) - 1

         if (listbc(3,ispot) .eq. 1) then
c           ::::: Cell i,j is on right side of a fine grid
            do 100 ivar = 1, nvar
               alloc(inlist + ivar) = -xfluxp(i,j,ivar)*dtc*hy
100         continue
         endif

         if (listbc(3,ispot) .eq. 2) then
c           ::::: Cell i,j on bottom side of fine grid
            do 200 ivar = 1, nvar
               alloc(inlist + ivar) = -yfluxm(i,j+1,ivar)*dtc*hx
200         continue
         endif

         if (listbc(3,ispot) .eq. 3) then
c           ::::: Cell i,j on left side of fine grid
            do 300 ivar = 1, nvar
               alloc(inlist + ivar) = -xfluxm(i+1,j,ivar)*dtc*hy
300         continue
         endif

         if (listbc(3,ispot) .eq. 4) then
c           ::::: Cell i,j on top side of fine grid
            do 400 ivar = 1, nvar
               alloc(inlist + ivar) = -yfluxp(i,j,ivar)*dtc*hx
400         continue
         endif
c
      ispot = ispot + 1
      if (ispot .gt. maxsp) go to 99
      go to 10
c
 99   return
      end
