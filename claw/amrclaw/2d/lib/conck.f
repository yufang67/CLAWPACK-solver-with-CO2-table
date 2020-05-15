c
c -----------------------------------------------------------
c
      subroutine conck(level, nvar)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      iadd(i,j,ivar)  = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iaddaux(i,j) = locaux + i - 1 + mitot*(j-1) +
     .                        mitot*mjtot*(mcapa-1)
c
c ******************************************************************
c conck - conservation check  for specified level
c         mostly a debugging tool
c         this assumes grids don't overlap
c ******************************************************************
c
c
c  grid loop for given level
c
      hx      = hxposs(level)
      hy      = hyposs(level)
      dt      = possk(level)
      totmass = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(store1,mptr)
         locaux = node(storeaux,mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost
c
         if (mcapa .eq. 0) then
           do 50 j  = nghost+1, mjtot-nghost
           do 50 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(i,j,1)) 
 50           continue
          else
c          # with capa array:
           do 60 j  = nghost+1, mjtot-nghost
           do 60 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(i,j,1))*alloc(iaddaux(i,j)) 
 60           continue
          endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    totmass = totmass * hx * hy
       write(outunit,*)" total mass = ", totmass
c
 99   return
      end
