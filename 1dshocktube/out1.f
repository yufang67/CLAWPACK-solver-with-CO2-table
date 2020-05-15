c
c
c =========================================================
      subroutine out1(maxmx,meqn,mbc,mx,xlower,dx,q,t,iframe,aux,maux)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 1 dimension
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>
c
c
      USE properties
      USE Interp_table
c      
      implicit double precision (a-h,o-z)
      dimension q  (1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, maux)
      character*10 fname1, fname2, fname3, fname4
      logical outaux, out_W

      common /comic/ rhol, Tl
      common /time/ time_initial
      
      
      CALL cpu_time(time)
      
      outaux = .false.
      out_W  = .true.
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.  
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file
c
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         fname3 = 'fort.axxxx'
         fname4 = 'SIN32.Cxxxx'
         nstp = iframe
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            fname3(ipos:ipos) = char(ichar('0') + idigit)
            fname4(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

         open(unit=50,file=fname1,status='unknown',form='formatted')
         open(unit=60,file=fname2,status='unknown',form='formatted')

c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:
      ngrids = 1
      mptr   = 1
      level  = 1

      write(50,1001) time-time_initial, mptr, level, mx
 1001 format(e20.10,'                CPU time',    /,
     &       i5,   '                 grid_number', /,
     &       i5,   '                 AMR_level',   /,
     &       i5,   '                 mx')

      write(50,1002) xlower,dx
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    dx', /)


        do 10 i=1,mx
          do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,m)) .lt. 1d-99) q(i,m) = 0.d0
             enddo
c
          write(50,1005)  (q(i,m), m=1,meqn)
 1005     format(6e18.8)
c
 10       continue
        write(50,*) ' '
 20     continue
      write(50,*) ' '

c
c ----------------------------------------------------------------------
c EDITED TO USE XMGRACE (03/11/2016) 
c
      if (out_W) then      	
c        # also output the aux arrays:
         open(unit=72,file=fname4,status='unknown',form='formatted')
         do 112 i=1,mx    
	     r    = q(i,1)
	     u    = q(i,2) / r
	     E    = q(i,3) / r	     
	     eint = E - 0.5d0*u*u - 5.0d6  
c
c            values from TABLE CO2
             IF (i.EQ.1)   CALL pressure(Tl,1d0/rhol,p_guess)
 
             CALL CO2BLLT_EQUI(p, T_v, c, x_v, a_v, res,
     &                         eint, 1d0/r, p_guess) 
             p_guess = p           
c    
             if (dabs(p) .lt. 1d-99) p = 0.d0           
c
            write(72,1006) p,r,c, x_v, T_v, eint, a_v, u, dx*(i-0.5d0) 
 1006     format(11e18.8)  
c
  112       continue
         write(72,*) ' '
         close(unit=72)
         endif
c
c ----------------------------------------------------------------------
c
      if (outaux) then 
c        # also output the aux arrays:
         open(unit=70,file=fname3,status='unknown',form='formatted')
         write(70,1001) mptr,level,mx
         write(70,1002) xlower,dx

c
         do 110 i=1,mx
            do m=1,maux
c              # exponents with more than 2 digits cause problems reading
c              # into matlab... reset tiny values to zero:
               if (dabs(aux(i,m)) .lt. 1d-99) aux(i,m) = 0.d0
            enddo
c
            write(70,1005) (aux(i,m), m=1,maux)
c
  110       continue
         write(70,*) ' '
         close(unit=70)
         endif

      write(60,1000) t,meqn,ngrids,maux

 1000 format(e26.16,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 maux'/,/)
c

      close(unit=50)
      close(unit=60)

      return
      end
