c
c
c =========================================================
      subroutine out1(maxmx,meqn,mbc,mx,xlower,dx,q,t0,iframe,aux,maux)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 1 dimension
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.
c     # Here its adapted to output just the single grid.
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>
c
c
      USE properties
      USE peng_robinson
      USE stiffened
      USE Interp_table
      USE solver_eos
      USE location, ONLY: phaseloca
c      
      implicit double precision (a-h,o-z)
      dimension q  (1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, maux)
      character*10 fname1, fname2, fname3, fname4
      logical outaux, out_W
      integer exitflag
c      common /comic/ rhol
      common /time/ time_initial
      common /e_cst/ e_const, Tguess
      common /flags/ i_flag,i_flag_phase
c      
c
c
      CALL cpu_time(time)
      
      outaux = .false.
      out_W  = .true.
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.  
c     # Here its adapted to output just the single grid.
c
c     # first create the file name and open file
c
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         fname3 = 'fort.axxxx'
         fname4 = 'HLLC.lxxxx'
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
	     eint = E - 0.5d0*u*u   
c
c   ------------------------------
c    EOS -->  p, c, x_v, T_v, a_v
c   ------------------------------
c
c          print* , 'flag_eos in out.f', i_flag

      IF (i_flag==1) THEN
c  TABLE CO2
          eint = eint - e_const   
c     print*, "CO2 TABLE IN out1.f"
             IF (i.EQ.1) pguess_out = 6.0e6    
c       CALL pressure(Tl,1d0/rhol,pguess_out)
c     print*, "CO2 TABLE interpolation Start",iframe,r,eint,pguess_out 
             CALL CO2BLLT_EQUI(p, T_v, c, x_v, a_v, res,
     &                         eint, 1d0/r, pguess_out) 
c     print*, "CO2 TABLE interpolation finishe",iframe,p,T_v    
               pguess_out = p
               i_flag_phase = res
c     print*, 'p_guess= ', press           
c     print*, 'x_v= ', x_v
c
c
c
      ELSEIF (i_flag==2) THEN
c     
        v = 1.0d0 / r
c     print*, eint
        energy = eint - 241.3d3
        CALL phaseloca(x_v, a_v, i_flag_phase, energy,v)
c      
c     print*, i_flag_phase,Tguess
      IF (i.EQ.1) Tguess=300.0
        CALL eos_1d(2, T_v, out_2, resnorm, Niter,
     &            exitflag, eint, Tguess, v, out3)
         Tguess = T_v
c     print*, T_v
        CALL pressure_pr(T_v,v,p)
        CALL soundspeed_pr(T_v,v,c)
c     print*, 'out1.f at i = ', i
c
c      
      ELSEIF (i_flag==3) THEN
c
         v = 1.0d0 / r
         energy = eint - e_const 
        CALL phaseloca(x_v, a_v, i_flag_phase, energy,v)
c
c      print*,i_flag_phase
        CALL pressure_st(i_flag_phase,v,energy, p)
        CALL temperature_st(i_flag_phase,v,energy,T_v)     
        CALL sound_speed_st(i_flag_phase,v,energy,c)
c
c
c
      ELSEIF (i_flag==4) THEN
c
         v = 1.0d0 / r
         energy = eint - e_const   
         CALL phaseloca(x_v, a_v, i_flag_phase, energy,v)
c     
        IF (i_flag_phase==5) THEN
           IF (i.EQ.1) pguess_out = 6.0e6
           CALL CO2BLLT_EQUI(p, T_v, c, x_v, a_v, res,
     &                     energy, v, pguess_out)
           pguess_out = p
           i_flag_phase = res
        ELSE
c         IF (i.EQ.1) Tguess=Tl
c         IF (i.EQ.1) print*,Tguess
           CALL eos_1d(4, T_v, out_2, resnorm, Niter,
     &            exitflag, energy, Tguess, v, out3)
           Tguess = T_v+1.0
           CALL pressure(T_v,v,p)
           CALL sound_speed(T_v,v,c)
        ENDIF
      ENDIF
c
c           
c
      if (dabs(p) .lt. 1d-99) p = 0.d0           
c
           write(72,1006) p,r,c, x_v, T_v, eint, a_v, u, dx*(i-0.5d0)
c     &     ,i-1.0d0+1.0d0, real(lag_phase) 
 1006     format(11e18.8)  
c
  112       continue
         write(72,*) ' '
         close(unit=72)
         endif
         CALL cpu_time(time_fin)
         print*,"SAVE HLLC FINISHIED", time_initial,time_fin
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

      write(60,1000) t0,meqn,ngrids,maux

 1000 format(e26.16,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 maux'/,/)
c

      close(unit=50)
      close(unit=60)

      return
      end
