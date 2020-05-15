c
c     =====================================================
      subroutine restart(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q)
c     =====================================================
c
c     # Set initial conditions for q.
      USE Grid
c      USE properties
c      USE peng_robinson, ONLY: interenergy_pr
c      USE stiffened
      USE var_const, ONLY: e_const, guessP,eguess_out,vguess_in
c    &                    ,Tguess,guessP_rpn
c
c
c
c
      implicit double precision (a-h,o-z)
c
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      character*10 fname1, fname2
      common /restrt_block/ tinitial, iframe

      CALL cpu_time(time_1)
      CALL MAKE_GRID()
      CALL cpu_time(time_2)



      iunit = 16

c     # first create the file name and open file
      fname1 = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
      open(iunit,file=fname1)

c     # Read grid parameters.
      read(iunit,*) igrid
      read(iunit,*) level
      read(iunit,*) mx_in
      read(iunit,*) my_in
      read(iunit,*) xlow_in
      read(iunit,*) ylow_in
      read(iunit,*) dx_in
      read(iunit,*) dy_in
      read(iunit,*)       
c     # Test for compatibility of grid resolution.
      if (mx_in .ne. mx .or. my_in .ne. my) then
         stop 'rstart.f : data not compatible'
      endif

c     # Read variables in from old fort.qXXXX file.
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            read(iunit,*) (q(i,j,m), m=1,meqn)
c            print*,i,j,'q', (q(i,j,m),m=1,meqn)
         enddo
         read(iunit,*)         
      enddo
!         read(iunit,*)
      close(iunit)
c  INITIAL GUESS 
      do 150 j=1-mbc,my+mbc
c          guessp_BC(j) = 0.6d6
          vguess_in(j)  = 1.0/500.0
          eguess_out(j) = -108.56d3 
          do 150 i=1-mbc,mx+mbc
c                 guessP_rpn(i,j) = 7d6
c                 Tguess(i,j) = 400
c                 c_BC(i,j)   = 194.0
c                 guessP(i,j) = 1.0d6
c                 c_BC(i,j)   = 239.0
  150        continue
c     # Read initial time in from fort.tXXXX file.      
      fname2 = 'fort.t'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
      open(iunit,file=fname2)
      read(iunit,*) tinitial
      close(iunit)


      write(*,*) 'Restarting from old output file ', fname1
      write(*,*) 'Simulation will be restarted at time t = ', tinitial
      write(*,*) 'Inital condition will not be output to a matlab ',
     &     'plot file'
      write(*,*) 

      return
      end
