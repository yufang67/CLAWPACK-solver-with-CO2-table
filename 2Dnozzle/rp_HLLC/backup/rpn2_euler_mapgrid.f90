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
    USE Interp_table, ONLY: CO2BLLT_EQUI
    USE var_const, ONLY: e_const,flag_diagno,guessP

    implicit none

    integer :: ixy, maxm, meqn, mwaves, mbc, mx, my, maux
    double precision ::    ql(1-mbc:maxm+mbc,meqn)
    double precision ::    qr(1-mbc:maxm+mbc,meqn)
    double precision ::     s(1-mbc:maxm+mbc,mwaves)
    double precision ::  wave(1-mbc:maxm+mbc,meqn,mwaves)
    double precision ::  amdq(1-mbc:maxm+mbc,meqn)
    double precision ::  apdq(1-mbc:maxm+mbc,meqn)
    double precision ::  auxl(1-mbc:maxm+mbc,*)
    double precision ::  auxr(1-mbc:maxm+mbc,*)

    double precision :: ql_state(4), qr_state(4)
    double precision :: rhol, ul, el, cl, pl, vl,Tl,xl,al,resl
    double precision :: rhor, ur, er, cr, pr, vr,Tr,xr,ar,resr
    double precision :: rhsqrtl, rhsqrtr, rhsq2
    double precision :: u, enth, delta(4), rho, v
    double precision :: wave_local(3,4), s_local(3),uv(2)
    double precision :: speeds(2,3), u2v2l, u2v2r,pguess,ccm
    integer :: i, m, mw, mu, mv, ixy1, info,jj
    logical :: efix

    integer :: mcapa,locrot, locarea
    double precision :: rot(4), area


!    data efix /.true./
    data efix /.false./

    maux   = 7
 

    call get_aux_locations_n(ixy,mcapa,locrot,locarea)

    do i = 2-mbc,mx+mbc-1
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
        IF (ixy==1) THEN
!
           CALL CO2BLLT_EQUI(pl, Tl, cl, xl, al, resl,&
     &                         el, 1.0/rhol, guessP(i-1,jj))
           CALL CO2BLLT_EQUI(pr, Tr, cr, xr, ar, resr,&
     &                         er, 1.0/rhor, guessP(i,jj))
           pguess = 0.5*(guessP(i-1,jj)+guessP(i,jj))
!
        ELSEIF (ixy==2) THEN
!
           CALL CO2BLLT_EQUI(pl, Tl, cl, xl, al, resl,&
     &                         el, 1.0/rhol, guessP(jj,i-1))
           CALL CO2BLLT_EQUI(pr, Tr, cr, xr, ar, resr,&
     &                         er, 1.0/rhor, guessP(jj,i))
           pguess = 0.5*(guessP(jj,i-1)+guessP(jj,i))
!
        ENDIF
!
!        # Get Roe averaged values
        rhsqrtl = sqrt(rhol)
        rhsqrtr = sqrt(rhor)
        rhsq2 = rhsqrtl + rhsqrtr
!
        uv(1) = (ul*rhsqrtl + ur*rhsqrtr) / rhsq2
        uv(2) = (vl*rhsqrtl + vr*rhsqrtr) / rhsq2
!### total enthalpy
        enth = ( (ql_state(4) + pl)/rhsqrtl + (qr_state(4) + pr)/rhsqrtr )/rhsq2
        ccm  = (cl*rhsqrtl + cr*rhsqrtr) / rhsq2
!  
!
        do m = 1,meqn
            delta(m) = qr_state(m) - ql_state(m)
        enddo

        ixy1 = 1
        call roe_solver(ixy1,uv,enth,delta,wave_local,s_local,info,ccm)

        if (info /= 0) then
            write(6,*) 'ixy = ', ixy
            write(6,*) 'enth = ', enth
            write(6,*) 'Called from rpn2 '
            write(6,*) ' '
            do m = 1,meqn
                write(6,'(2E24.16)') ql_state(m), qr_state(m)
            enddo
            write(6,*) ' '
            stop
        endif

        do mw = 1,mwaves
            speeds(1,mw) = min(s_local(mw),0.d0)
            speeds(2,mw) = max(s_local(mw),0.d0)
        enddo

        if (efix) then
!           # This modifies the speeds, but we will still have
!           # s(mw) = speeds(mw,1) + speeds(mw,2)
            call apply_entropy_fix(ql_state,qr_state,cl,cr, &
            wave_local, speeds, pguess)
        endif

        area = auxl(i,locarea)
        do mw = 1,mwaves
            call rotate2_tr(rot,wave_local(mw,2),wave_local(mw,3))
            speeds(1,mw) = area*speeds(1,mw)
            speeds(2,mw) = area*speeds(2,mw)
        enddo

        do m = 1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw = 1,mwaves
                wave(i,m,mw) = wave_local(mw,m)
                s(i,mw) = speeds(1,mw) + speeds(2,mw)
                amdq(i,m) = amdq(i,m) + speeds(1,mw)*wave(i,m,mw)
                apdq(i,m) = apdq(i,m) + speeds(2,mw)*wave(i,m,mw)
            enddo
        enddo
    enddo


    return
    end subroutine rpn2

    subroutine apply_entropy_fix(ql,qr,cl, cr,wave_local,speeds,pguess)

    USE Interp_table, ONLY: CO2BLLT_EQUI
    USE var_const, ONLY: e_const,flag_diagno
    implicit none

    double precision :: ql(4), qr(4),wave_local(3,4)
    double precision :: speeds(2,3), s1, s2, s3, s4
    double precision :: sl, sml, smr, sr, qml(4),qmr(4)
    double precision :: ul, cl, ur, cr, pml, pmr, cml,cmr
    double precision :: sfract
    double precision :: press, temp,sound,qual,mass,res,rhoml,eml,uml,vml,pguess
    double precision :: rhomr,emr,umr,vmr
    logical :: trans1, trans3

    integer :: m


    s1 = speeds(1,1) + speeds(2,1)
    s2 = speeds(1,2) + speeds(2,2)
    s3 = speeds(1,3) + speeds(2,3)

    do m = 1,4
        qml(m) = ql(m) + wave_local(1,m)
    enddo
    sl = ql(2)/ql(1) - cl
    rhoml = qml(1)
    uml   = qml(2)/qml(1)
    vml   = qml(3)/qml(1)
    eml   = qml(4)/qml(1) - 0.5*(uml*uml+vml*vml) - e_const

    CALL CO2BLLT_EQUI(press, temp, sound, qual, mass, res,&
  &                   eml, 1.0/rhoml, pguess)

    pml = press
    sml = qml(2)/qml(1) - sound

!     # Check the 1-wave
    trans1 = .false.
    if (sl < 0 .AND. sml > 0) then
    !        # apply transonic entropy fix
        trans1 = .true.
        sfract = (sml - s1)/(sml - sl)
        speeds(1,1) = sfract*sl
        speeds(2,1) = (1-sfract)*sml
    endif

!     # If the 1-wave is transonic,then we are done...
!     # Otherwise, we have to check the 3-wave
    if ( .NOT. trans1) then
        do m = 1,4
            qmr(m) = qr(m) - wave_local(3,m)
        enddo
    rhomr = qmr(1)
    umr   = qmr(2)/qmr(1)
    vmr   = qmr(3)/qmr(1)
    emr   = qmr(4)/qmr(1) - 0.5*(umr*umr+vmr*vmr) - e_const
    sr    = qr(2)/qr(1) + cr
    CALL CO2BLLT_EQUI(press, temp, sound, qual, mass, res,&
  &                   emr, 1.0/rhomr, pguess)


    pmr = press
    smr = qmr(2)/qmr(1) + sound
    trans3 = .false.
        if (smr < 0 .AND. sr > 0) then
!           # apply transonic entropy fix
            trans3 = .true.
            sfract = (sr - s3)/(sr - smr)
            speeds(1,3) = sfract*smr
            speeds(2,3) = (1-sfract)*sr
         endif
    endif

!    if (trans1) then
!        write(6,*) '1-wave is transonic'
!    elseif (trans3) then
!        write(6,*) '3-wave is transonic'
!    endif

    end subroutine apply_entropy_fix
