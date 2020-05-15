subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,&
&               ql,qr,auxl,auxr,k,          &
&               wave,s,amdq,apdq)
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
!                                    and right state ql(i,:) 
! Solve Riemann problems along one slice of data.
!      USE Interp_table
!!      USE peng_robinson
!      USE solver_eos
!!      USE location
!!      USE stiffened
!      USE properties
!      USE var_const, ONLY: Tguess,e_const,guessP, flag_diagno,&
!&                          c_BC,i_flag
!
    implicit none
!
    integer,intent(in) :: ixy, maxm, meqn, mwaves, mbc, mx, k
    double precision,intent(in)  ::    ql(1-mbc:maxm+mbc,meqn)
    double precision,intent(in)  ::    qr(1-mbc:maxm+mbc,meqn)
    double precision,intent(out) ::     s(1-mbc:maxm+mbc,mwaves)
    double precision,intent(out) ::  wave(1-mbc:maxm+mbc,meqn,mwaves)
    double precision,intent(out) ::  amdq(1-mbc:maxm+mbc,meqn)
    double precision,intent(out) ::  apdq(1-mbc:maxm+mbc,meqn)
    double precision,intent(in)  ::  auxl(1-mbc:maxm+mbc,*)
    double precision,intent(in)  ::  auxr(1-mbc:maxm+mbc,*)
!
    double precision :: ql_state(4), qr_state(4)
    double precision :: rhol, ul, el, cl, pl, vl
    double precision :: rhor, ur, er, cr, pr, vr
    double precision :: qml(4), qmr(4)
    real(8)::  press, sound,xr, xl, a_out, res,Tl, Tr, ar, al
    real(8)::  guess_t,res1,res2!,pguess
!
    real(8):: sr,sl,sm


    double precision ::  u2v2l, u2v2r
    integer :: i, m, mw, ixy1, info

    integer :: locarea
    double precision :: rot(4), area
!    
    DO i = 1,mx
!Construct the rotation matrix R
        IF(ixy==1) THEN
          rot(1) = auxl(i,1)
          rot(2) = auxl(i,2)
          rot(3) = -rot(2)
          rot(4) =  rot(1)
          locarea = 3
        ELSEIF (ixy==2) THEN
          rot(1) = auxl(i,4)
          rot(2) = auxl(i,5)
          rot(3) = rot(2)  !??
          rot(4) = -rot(1) !??
          locarea = 6
        ENDIF

!        print*, rot(1),rot(2),rot(3),rot(4)
       
!Compute rotational q for computational space(only u or v )
        do m = 1,meqn
            qr_state(m) = ql(i,m)
            ql_state(m) = qr(i-1,m)
        enddo
        call rotate2(rot,ql_state(2),ql_state(3))
        call rotate2(rot,qr_state(2),qr_state(3))
!        
        rhol = ql_state(1)
        rhor = qr_state(1)
!        
        ul = ql_state(2)/rhol
        ur = qr_state(2)/rhor
!
        vl = ql_state(3)/rhol
        vr = qr_state(3)/rhor
! 
        u2v2l = ul*ul + vl*vl
        u2v2r = ur*ur + vr*vr
!
!
        el = ql_state(4)/rhol - 0.5*u2v2l      !internal energy
        er = qr_state(4)/rhor - 0.5*u2v2r
!        print*,i,'left=', rhol,ul,vl,el
!        print*,i,k,'right=',rhor,ur,vr,er
!-------------------------------------------------------------------------
!  EOS:  Using er,el, rhor, rhol to find press_r and press_l,c_r,c_l
!      print*,'#values from TABLE CO2# in rpn',i
!         
          pl = (1.4-1.0)*el*rhol
          cl = sqrt(pl / rhol * 1.4)
!
!
          pr = (1.4-1.0)*er*rhor
          cr = sqrt(pr / rhor * 1.4)
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
         DO mw=1,mwaves
            IF (s(i,mw)/=s(i,mw)) THEN
               Print*, 'HLLC solver ',' i=',i, 'mw=',mw, 's=',s(i,mw)
               Print*, 'ul=',ul,'ur=',ur,'cl=',cl,'cr=',cr
               STOP
            ENDIF
         ENDDO
!
!        intermediate states
!
         qml(1) = 1.0
         qml(2) = sm
         qml(3) = vl
         qml(4) = (ql_state(4)/rhol + (sm - ul)*(sm + pl/   &
     &                        (rhol*(sl - ul))))
  
!         qml(1:4) = qml(1:4) * rhol*(sl - ul)/(sl - sm)
!
         qmr(1) = 1.0
         qmr(2) = sm
         qmr(3) = vr
         qmr(4) = (qr_state(4)/rhor + (sm - ur)*(sm + pr/     &
     &                        (rhor*(sr - ur))))
!         qmr(1:4) = qmr(1:4) * rhor*(sr - ur)/(sr - sm)
!
!        # Compute the waves.
!
         DO m =1,4
            qmr(m) = qmr(m) * rhor*(sr - ur)/(sr - sm)
            qml(m) = qml(m) * rhol*(sl - ul)/(sl - sm)
            wave(i,m,1) = qml(m)  - ql_state(m)
            wave(i,m,2) = qmr(m)  - qml(m)
            wave(i,m,3) = qr_state(m) - qmr(m)
            IF (wave(i,m,1)/=wave(i,m,1)) THEN
               Print*, 'HLLC solver ','i ', i,'m ',m,'waves ',wave(i,m,1)
               Print*, 'qml',k,ixy, pl,pr
               STOP
             ENDIF

         ENDDO
    ENDDO
!
!     # compute Godunov flux f0:
!     --------------------------
!
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves
!
      DO  m=1,4
         DO  i=1-mbc, mx+mbc
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
! CHECK
       DO i = 1-mbc,mx+mbc
           DO m = 1,4
             IF ( (amdq(i,m)/=amdq(i,m)) .OR. (apdq(i,m)/=apdq(i,m)) ) THEN
               Print*, 'HLLC solver ',i,'amdq', amdq(i,m),'apdq', apdq(i,m)
               STOP
             ENDIF
           ENDDO
       ENDDO
    return
    end subroutine rpn2
