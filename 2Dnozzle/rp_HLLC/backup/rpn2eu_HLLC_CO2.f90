subroutine rpn2(ixy,maxm, meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! HLLC-solver for the Euler equations with mapped grids

! waves: 3
! equations: 4

! Conserved quantities:
!       1 density
!       2 x_momentum
!       3 y_momentum
!       4 energy
!
! Auxiliary quantities:
!       1 nx(i-0.5,j)
!       2 ny(i-0.5,j)
!       3 length_ratio_left ratio of length of left face to dyc
!       4 nx(i,j-0.5)
!       5 ny(i,j-0.5)
!       6 length_ratio_bottom   ratio of length of bottom face to dxc
!       7 cell_area   ratio of cell area to dxc*dyc
! Rotate Ql Qr in physical space to computational space
!       Ql*R , Qr*R
!
! Compute p,c by using e and rho
!
! Compute s, wave, intermidate states qml, qmr
!
! Compute modified s 
!         rotational wave by wave * R^-1           
!         modified amdq and apdq
!
!
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:) ??
! Solve Riemann problems along one slice of data.
      USE Interp_table
      USE peng_robinson
      USE solver_eos
      USE location
      USE stiffened
      USE properties
!
    implicit none
!
    integer,intent(in) :: ixy, maxm, meqn, mwaves, mbc, mx
    double precision,intent(in)  ::    ql(1-mbc:maxm+mbc,meqn)
    double precision,intent(in)  ::    qr(1-mbc:maxm+mbc,meqn)
    double precision,intent(out) ::     s(1-mbc:maxm+mbc,mwaves)
    double precision,intent(out) ::  wave(1-mbc:maxm+mbc,meqn,mwaves)
    double precision,intent(out) ::  amdq(1-mbc:maxm+mbc,meqn)
    double precision,intent(out) ::  apdq(1-mbc:maxm+mbc,meqn)
    double precision,intent(in)  ::  auxl(1-mbc:maxm+mbc,maux)
    double precision,intent(in)  ::  auxr(1-mbc:maxm+mbc,maux)
!
    real(8),parameter :: e_const=0.5d6
    real(8),parameter :: Tguess = 300.0
    integer,parameter :: i_flag = 1 
!
    double precision :: ql_state(4), qr_state(4)
    double precision :: rhol, ul, el, cl, pl, vl
    double precision :: rhor, ur, er, cr, pr, vr
    double precision :: qml(4), qmr(4)
    integer exitflag
    real(8)::  press, sound,xr, xl, a_out, res,Tl, Tr, ar, al
    real(8)::  guess_t,pguess,res1,res2
!
    real(8):: sr,sl,sm


    double precision ::  u2v2l, u2v2r
    integer :: i, m, mw, mu, mv, ixy1, info

    integer :: mcapa,locrot, locarea,maux
    double precision :: rot(4), area

    maux = 7


    call get_aux_locations_n(ixy,mcapa,locrot,locarea)
!     print*, ixy,mcapa,locrot,locarea
    DO i = 2-mbc,mx+mbc
!Construct the rotation matrix R
        rot(1) = auxl(i,locrot)
        rot(2) = auxl(i,locrot+1)
        call compute_tangent(rot)
!        print*, rot(1),rot(2),rot(3),rot(4)
!Compute rotational q for computational space(only u or v )
        do m = 1,meqn
            qr_state(m) = qr(i,m)
            ql_state(m) = ql(i-1,m)
        enddo
        call rotate2(rot,ql_state(2),ql_state(3))
        call rotate2(rot,qr_state(2),qr_state(3))

        rhol = ql_state(1)
        rhor = qr_state(1)

        ul = ql_state(2)/rhol
        ur = qr_state(2)/rhor

        vl = ql_state(3)/rhol
        vr = qr_state(3)/rhor
        
        u2v2l = ul*ul + vl*vl
        u2v2r = ur*ur + vr*vr

        el = ql_state(4)/rhol - 0.5*u2v2l     !internal energy
        er = qr_state(4)/rhor - 0.5*u2v2r
!        print*,i,'left=', rhol,ul,vl,el
!        print*,i,'right=',rhor,ur,vr,er
!-------------------------------------------------------------------------
!  EOS:  Using er,el, rhor, rhol to find press_r and press_l,c_r,c_l
      IF (i_flag==1)THEN                           ! #values from TABLE CO2#
         pguess = 2.0d6
         el = el - e_const
         
         CALL CO2BLLT_EQUI(press, Tl, sound, xl, al, res1,&
     &                 el,1.d0/rhol,pguess)
!        
         pl = press
         cl = sound
         pguess = press
!         print*, 'left state finiched in rpn2'
         er = er - e_const
         CALL CO2BLLT_EQUI(press, Tr, sound, xr, ar, res2,&
     &                 er,1.d0/rhor,pguess)
         pr = press
         cr = sound
         pguess = press
!         print*, 'right state finiched in rpn2'
!         Print*, 'NODE', i 
!!        IF (i==500 .OR. i==501) THEN
!!     print*,i, res1, res2
!!        ENDIF

!      ELSE IF (i_flag==2) THEN
!         energy = el - 241.3d3
!         CALL phaseloca(xl, al, i_flag_phase, energy,vl)
!         IF (i_flag_phase==1 .OR. i_flag_phase==5) THEN
!            PRINT*, 'liquid presents in the left wave in rp1.f at',i
!!            STOP
!         ENDIF
!
!         IF (i.EQ.0) guess_t=300.0+1.0
!         CALL eos_1d(2, Tl, out_2, resnorm, Niter,
!     &            exitflag, el, guess_t, vl, out3)
!         guess_t = Tl+1.0
!         CALL pressure_pr(Tl,vl,press)
!         CALL soundspeed_pr(Tl,vl,sound)
!         press_l = press
!         cl = sound
!!     print*, "i node left",i, Tl,press_l,cl
!         energy = er - 241.3d3
!         CALL phaseloca(xr, ar, i_flag_phase, energy,vr)
!         IF (i_flag_phase==1 .OR. i_flag_phase==5) THEN
!            PRINT*, 'liquid presents in the right wave in rp1.f at',i
!!            STOP
!         ENDIF
!!
!         CALL eos_1d(2, Tr, out_2, resnorm, Niter,
!     &            exitflag, er, guess_t, vr, out3)
!         guess_t = Tr-1.0
!         CALL pressure_pr(Tr,vr,press)
!         CALL soundspeed_pr(Tr,vr,sound)
!         press_r = press
!         cr = sound
!!     print*, 'temperature', Tl, Tr, guess_t
!!     print*, "i node right",i, Tr,press_r,cr
!     ELSE IF (i_flag==3) THEN
!         energy = el - e_const
!         CALL phaseloca(xl, al, i_flag_phase, energy,vl)
!!
!!      print*, 'flag_phase=', i_flag_phase,i
!         CALL pressure_st(i_flag_phase,vl,energy, press_l)
!         CALL sound_speed_st(i_flag_phase,vl,energy,cl)
!!
!         energy = er - e_const
!         CALL phaseloca(xr, ar, i_flag_phase, energy,vr)
!!
!!      print*, 'flag_phase=', i_flag_phase,i
!        CALL pressure_st(i_flag_phase,vr,energy, press_r)
!         CALL sound_speed_st(i_flag_phase,vr,energy,cr)
!!
!!
!      ELSE IF (i_flag==4) THEN
!         energy = el - e_const
!         CALL phaseloca(xl, al, i_flag_phase, energy,vl)
!!     
!        IF (i_flag_phase==5) THEN
!           IF (i.EQ.1) pguess_out = 6.0e6
!           CALL CO2BLLT_EQUI(press_l, Tl, cl, xl, al, res1,
!     &                     energy, vl, pguess_out)
!           pguess_out = press_l
!           i_flag_phase = res1
!        ELSE
!!         IF (i.EQ.1) Tguess=Tl
!           CALL eos_1d(4, Tl, out_2, resnorm, Niter,
!     &            exitflag, energy, Tguess, vl, out3)
!           Tguess = Tl
!           CALL pressure(Tl,vl,press_l)
!           CALL sound_speed(Tl,vl,cl)
! ENDIF
!!
!        energy = er - e_const
!         CALL phaseloca(xr, ar, i_flag_phase, energy,vr)
!!     
!        IF (i_flag_phase==5) THEN
!           IF (i.EQ.1) pguess_out = 6.0e6
!           CALL CO2BLLT_EQUI(press_r, Tr, cr, xr, ar, res2,
!     &                     energy, vr, pguess_out)
!           pguess_out = press_r
!           i_flag_phase = res2
!        ELSE
!!         IF (i.EQ.1) Tguess=Tl
!          CALL eos_1d(4, Tr, out_2, resnorm, Niter,
!    &            exitflag, energy, Tguess, vr, out3)
!          Tguess = Tr
!           CALL pressure(Tr,vr,press_r)
!           CALL sound_speed(Tr,vr,cr)
!        ENDIF
!
      ENDIF
!      print*,i, 'interpolation in CO2 table finished'
!! 
! 
!        speed of sound  for signal velocity calculation
! 
         sl     = min(ul-cl, ur-cr)
         sr     = max(ul+cl, ur+cr)
         sm     = (pr-pl + rhol*ul*(sl- ul) - rhor*ur*(sr- ur))&
     &                  / (rhol*   (sl- ul) - rhor*   (sr- ur))
         area = auxl(i,locarea)
!         print*,area,auxl(i,3)
         s(i,1) = sl*area
         s(i,2) = sm*area
         s(i,3) = sr*area
!
!        intermediate states
!
         qml(1) = 1d0
         qml(2) = sm
         qml(3) = 1d0   !????????????????????
         qml(4) = (ql_state(4)/rhol + (sm - ul)*(sm + pl/   &
     &                        (rhol*(sl - ul))))
         qml(1:4) = qml(1:4) * rhol*(sl - ul)/(sl - sm)
!
         qmr(1) = 1d0
         qmr(2) = sm
         qmr(3) = 1d0   !???????????????????
         qmr(4) = (qr_state(4)/rhor + (sm - ur)*(sm + pr/     &
     &                        (rhor*(sr - ur))))
         qmr(1:4) = qmr(1:4) * rhor*(sr - ur)/(sr - sm)
!         print*,qml(1),qmr(1),qml(2),qmr(2)
!
!        # Compute the waves.
!
         wave(i,1:4,1) = qml(1:4)  - ql_state(1:4)
         wave(i,1:4,2) = qmr(1:4)  - qml(1:4)
         wave(i,1:4,3) = qr_state(1:4) - qmr(1:4)
!         print*, wave(i,:,:)
    ENDDO
!         print*,'s and wave completed'
!
!     # compute Godunov flux f0:
!     --------------------------
!
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves
!
      DO  m=1,4
         DO  i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            DO  mw=1,mwaves
            CALL rotate2_tr(rot,wave(i,2,mw),wave(i,3,mw))
!            print*, wave(i,2,mw),wave(i,3,mw)
               if (s(i,mw) .lt. 0.d0) then
                   amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
!                   print*,'amdq', amdq(i,m)
                 else
!                   print*,'apdq', apdq(i,m)
                   apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                 endif
             ENDDO
         ENDDO
       ENDDO
!
!   print*, 'rpn2 has been called'
    return
    end subroutine rpn2
