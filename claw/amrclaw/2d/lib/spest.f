c
c -------------------------------------------------------------
c
      subroutine spest (nvar,naux,lcheck,iflags,isize,jsize)
c
      implicit double precision (a-h,o-z)

      dimension  iflags (0:isize+1,0:jsize+1)
      include  "call.i"

 
c :::::::::::::::::::::::::: SPEST :::::::::::::::::::::::::::::::::::
c For all grids at level lcheck:
c   Call user-supplied routine flag2refine to flag any points where
c   refinement is desired based on user's criterion.  
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c   initialize iflags here, so can put user flags right into it.
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
          dx     = hxposs(lcheck)
          dy     = hyposs(lcheck)
          xleft  = rnode(cornxlo,mptr)
          ybot   = rnode(cornylo,mptr)
          xlow   = xleft - nghost*dx
          ylow   = ybot - nghost*dy
c
          locbig = igetsp(mitot*mjtot*nvar)
          node(tempptr,mptr) = locbig
c         # straight copy into scratch array so don't mess up latest soln.
          do 10 i = 1, mitot*mjtot*nvar
 10          alloc(locbig+i-1) = alloc(locnew+i-1)

          call bound(time,nvar,nghost,alloc(locbig),mitot,mjtot,mptr,
     1               alloc(locaux),naux)
c
c get user flags for refinement, which might be based on spatial gradient, 
c for example.  Use old values of soln at time t  before
c integration to get accurate boundary gradients
c
      if (tolsp .gt. 0.) then
         locamrflags = igetsp(mitot*mjtot)
	    do 20 i = 1, mitot*mjtot
 20         alloc(locamrflags+i-1) = goodpt

c        # call user-supplied routine to flag any points where 
c        # refinement is desired based on user's criterion.  
c        # Default version compares spatial gradient to tolsp.

         call flag2refine(nx,ny,nghost,nvar,naux,xlow,ylow,dx,dy,
     &              time,lcheck,tolsp,alloc(locbig),alloc(locaux),
     &              alloc(locamrflags), goodpt, badpt )

c        Put flags in iflags array now, so can reclaim space.
c        Note change of dimension of amrflags array:
c        3rd dim = 1 here, elsewhere flag array has idim3 = nvar
c
         idim3 = 1   ! 3rd dim = 1 here, elsewhere is nvar
         call setflags (iflags,isize,jsize,
     1                  alloc(locamrflags),idim3,mitot,mjtot,mptr)
         call reclam(locamrflags, mitot*mjtot)
      endif

c   previously used to reclam locbig space here. now save to reuse in errest
c   reclam locbig space afterwards.
c     call reclam(locbig,mitot*mjtot*nvar)

      mptr = node(levelptr,mptr)
      if (mptr .ne. 0) go to 5
c
      return
      end
