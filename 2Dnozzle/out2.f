c
c
c =========================================================
      subroutine out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,t,iframe,aux,maux)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>
c
c
c      USE solver_eos
c      USE location
c      USE properties
      USE Interp_table, ONLY: CO2BLLT_EQUI
      USE var_const, ONLY: guessP,e_const, flag_diagno
c      
      implicit double precision (a-h,o-z)
      dimension   q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      character*10 fname1, fname2, fname3,fname4
      logical outaux, out_W

      outaux = .false.
      out_W  = .true.
c
c     # first create the file name and open file
c
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         fname3 = 'fort.axxxx'
         fname4 = 'NOZ1.Cxxxx'
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
      mptr = 1
      level = 1

      write(50,1001) mptr,level,mx,my
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')

      write(50,1002) xlower,ylower,dx,dy
 1002 format(e26.16,'    xlow', /,
     &       e26.16,'    ylow', /,
     &       e26.16,'    dx', /,
     &       e26.16,'    dy',/)
c
      do 20 j=1-mbc,my+mbc
        do 10 i=1-mbc,mx+mbc
             do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,j,m)) .lt. 1d-99) q(i,j,m) = 0.d0
             enddo
c
          write(50,1005) (q(i,j,m), m=1,meqn)
 1005     format(4e26.16)
c
 10       continue
        write(50,*) ' '
 20     continue
      write(50,*) ' '
c============================================================================
c           08/2017
c===========================================================================
      IF (out_W) THEN
         OPEN(UNIT=72,FILE=fname4,STATUS='unknown',FORM='formatted')
         DO j=1-mbc,my+mbc
            DO  i=1-mbc,mx+mbc
                r   =  q(i,j,1)
                u   =  q(i,j,2)/r
                vv   =  q(i,j,3)/r
                E   =  q(i,j,4)/r
                eint=  E - 0.5d0*u*u - 0.5d0*vv*vv - e_const
c
c
               CALL CO2BLLT_EQUI(p, T_v, c, x_v, a_v, res,
     &                           eint, 1d0/r, guessP(i,j))
c
               IF (dabs(p) .lt. 1d-99) p = 0.d0
               IF (dabs(r) .lt. 1d-99) r = 0.d0
               IF (dabs(c) .lt. 1d-99) c = 0.d0
               IF (dabs(res) .lt. 1d-99) res = 0.d0
               IF (dabs(x_v) .lt. 1d-99) x_v = 0.d0
               IF (dabs(a_v) .lt. 1d-99) a_v = 0.d0
               IF (dabs(u) .lt. 1d-99) u = 0.d0
               IF (dabs(vv) .lt. 1d-99) vv = 0.d0
c
               xcposit = (i-0.5)*dx
               ycposit = (j-0.5)*dy
               CALL mapc2p(xcposit,ycposit,xp,yp)

          write(72,1006) p,r,c, x_v, T_v, eint, u,vv, 
     &                   xp,yp,a_v,res 
 1006     format(12e18.8)
c
            ENDDO
          write(72,*) ' '
         ENDDO
          write(72,*) ' '
          close(unit=72)
c
      ENDIF
c
c         print*,"SAVE solut FINISHIED"
c
c===========================================================================
c
      if (outaux) then 
c     # also output the aux arrays:
      open(unit=70,file=fname3,status='unknown',form='formatted')
      write(70,1001) mptr,level,mx,my
      write(70,1002) xlower,ylower,dx,dy
      do 120 j=1,my
         do 110 i=1,mx
            do m=1,maux
c              # exponents with more than 2 digits cause problems reading
c              # into matlab... reset tiny values to zero:
               if (dabs(aux(i,j,m)) .lt. 1d-99) aux(i,j,m) = 0.d0
            enddo
c
            write(70,1005) (aux(i,j,m), m=1,maux)
c
  110       continue
         write(70,*) ' '
  120    continue
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
