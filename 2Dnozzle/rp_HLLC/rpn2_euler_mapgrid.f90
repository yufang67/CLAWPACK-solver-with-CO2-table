subroutine rpn2(ixy,maxm, meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,jj,wave,s,amdq,apdq)
! Roe-solver for the Euler equations with mapped grids

! waves: 3
! equations: 4

! Conserved quantities:
!       1 density
!       2 x_momentum
!       3 y_momentum
!       4 energy

! Need to fill in the auxiliary variables

! Solve Riemann problems along one slice of data.
    USE def_constants, ONLY: pr
    USE Interp_table, ONLY: CO2BLLT_EQUI
    USE var_const, ONLY: e_const,flag_diagno,flag_out
    implicit none

    integer :: ixy, maxm, meqn, mwaves, mbc, mx, my, maux
    REAL(pr) ::    ql(1-mbc:maxm+mbc,meqn)
    REAL(pr) ::    qr(1-mbc:maxm+mbc,meqn)
    REAL(pr) ::     s(1-mbc:maxm+mbc,mwaves)
    REAL(pr) ::  wave(1-mbc:maxm+mbc,meqn,mwaves)
    REAL(pr) ::  amdq(1-mbc:maxm+mbc,meqn)
    REAL(pr) ::  apdq(1-mbc:maxm+mbc,meqn)
    REAL(pr) ::  auxl(1-mbc:maxm+mbc,*)
    REAL(pr) ::  auxr(1-mbc:maxm+mbc,*)

    REAL(pr) :: ql_state(4), qr_state(4)
    REAL(pr) :: rhol, ul, el, cl, pl, vl,Tl,xl,al,resl
    REAL(pr) :: rhor, ur, er, cr, prr, vr,Tr,xr,ar,resr

    REAL(pr) :: u, enth, delta(4), rho, v
    REAL(pr) :: wave_local(3,4), s_local(3),uv(2)
    REAL(pr) :: speeds(2,3), u2v2l, u2v2r,pguess
    integer :: i, m, mw, mu, mv, ixy1, info,jj    
    integer :: mcapa,locrot, locarea
    REAL(pr) :: rot(4), area

!
!
!    print*,'rpn start'
    maux   = 7
    flag_diagno = 2
 
    call get_aux_locations_n(ixy,mcapa,locrot,locarea)

    do i = 2-mbc,mx+2

        rot(1) = auxl(i,locrot)
        rot(2) = auxl(i,locrot+1)
        call compute_tangent(rot)

        do m = 1,meqn
            ql_state(m) = qr(i-1,m)
            qr_state(m) = ql(i,m)
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

        el = ql_state(4)/rhol - e_const - 0.5*u2v2l
        er = qr_state(4)/rhor - e_const - 0.5*u2v2r
!
!        IF (ixy==1) THEN
!
!        print*, 'rpn i=',i
!        print*, 'left in', el,1.0/rhol
           CALL CO2BLLT_EQUI(pl, Tl, cl, xl, al, resl,&
     &                         el, 1.0/rhol, 0.0d0)
!
!
!        print*, 'left out',pl,Tl,cl,xl,al,'flag',flag_out
!
!        print*, 'right in', er,1.0/rhor

           CALL CO2BLLT_EQUI(prr, Tr, cr, xr, ar, resr,&
     &                         er, 1.0/rhor, 0.0d0)
!        print*, 'right out', prr,Tr,cr,xr,ar,'flag',flag_out
!           pguess = 0.5*(guessP(i-1,jj)+guessP(i,jj))
!
!        ELSEIF (ixy==2) THEN
!!
!           CALL CO2BLLT_EQUI(pl, Tl, cl, xl, al, resl,&
!     &                         el, 1.0/rhol, guessP(jj,i-1))
!
!           CALL CO2BLLT_EQUI(pr, Tr, cr, xr, ar, resr,&
!     &                         er, 1.0/rhor, guessP(jj,i))
!           pguess = 0.5*(guessP(jj,i-1)+guessP(jj,i))
!
!        ENDIF
!
        IF ( (cl /= cl) .OR. (cr /= cr) ) THEN
            write(6,*) 'ixy = ', ixy,i
            write(6,*) 'cl = ', cl,'cr = ',cr
            write(6,*) 'Called from rpn2, sound speed'
            write(6,*) ' '
            write(6,*) 'left' , el, 1.0/rhol
            write(6,*) 'right', er, 1.0/rhor
!            do m = 1,meqn
!                write(6,'(2E24.16)') ql_state(m), qr_state(m)
!            enddo
            write(6,*) ' '
            stop
         ENDIF
!
!         print*, 'before hllc solver'
         CALL hllc_solver(pl,prr,ul,ur,vl,vr,cl,cr,rhol,rhor,&
&                          ql_state,qr_state,wave_local,s_local,info)
! 
!        print*, wave_local
        area = auxl(i,locarea)
        do mw = 1,mwaves
            call rotate2_tr(rot,wave_local(mw,2),wave_local(mw,3))
            s(i,mw) = area*s_local(mw)
        enddo

        do m = 1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw = 1,mwaves
                wave(i,m,mw) = wave_local(mw,m)
                if (s(i,mw) .lt. 0.d0) then
                   amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
!                   print*,'amdq', amdq(i,m)
                 else
!                   print*,'apdq', apdq(i,m)
                   apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                 endif
            enddo
        enddo
    enddo
!
!           print*, 'rpn finish'
    end subroutine rpn2
