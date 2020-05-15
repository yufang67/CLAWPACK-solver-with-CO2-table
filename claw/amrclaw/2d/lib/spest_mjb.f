c
c -------------------------------------------------------------
c
      subroutine spest (nvar,naux,lcheck,iflags,isize,jsize)
c
      implicit double precision (a-h,o-z)

      dimension  iflags (0:isize+1,0:jsize+1)
      include  "call.i"

 
c :::::::::::::::::::::::::: SPEST :::::::::::::::::::::::::::::::::::
c for all grids at level lcheck:
c  spatial error estimation (simple way for user to control
c  flagging, without Richardson.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c   initialize iflags here, so can put spatial error right into it.
c
       do 4 i = 1, isize
       do 4 j = 1, jsize
 4        iflags(i,j) = 0
c
       mptr = lstart(lcheck)
 5     continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locnew = node(store1,mptr)
          locaux = node(storeaux,mptr)
          time   = rnode(timemult,mptr)
c
          locbig = igetsp(mitot*mjtot*nvar)
          node(tempptr,mptr) = locbig
c         # straight copy into scratch array so don't mess up latest soln.
          do 10 i = 1, mitot*mjtot*nvar
 10          alloc(locbig+i-1) = alloc(locnew+i-1)

          call bound(time,nvar,nghost,alloc(locbig),mitot,mjtot,mptr,
     1               alloc(locaux),naux)
c
c get user estimate for refinement, could be spatial error, or moving
c allowable patch.  use old values of soln at time t  before
c integration to get accurate boundary gradients
c
      if (tolsp .gt. 0.) then
         locerrsp = igetsp(mitot*mjtot)
	 do 20 i = 1, mitot*mjtot
 20         alloc(locerrsp+i-1) = goodpt
         call errsp(alloc(locbig), alloc(locerrsp), mitot, mjtot,
     &              nvar, nghost, tolsp, time, goodpt, badpt )
c
c        put flags in iflags array now, so can reclaim space
c        note change of dimension than flag array from errest
c        3rd dim = 1 here, elsewhere is nvar
c
         idim3 = 1   ! 3rd dim = 1 here, elsewhere is nvar
         call setflags (iflags,isize,jsize,
     1                  alloc(locerrsp),idim3,mitot,mjtot,mptr)
         call reclam(locerrsp, mitot*mjtot)
      endif

c   previously used to reclam locbig space here. now save to reuse in errest
c   reclam locbig space afterwards.
c     call reclam(locbig,mitot*mjtot*nvar)

      mptr = node(levelptr,mptr)
      if (mptr .ne. 0) go to 5
c
      return
      end
