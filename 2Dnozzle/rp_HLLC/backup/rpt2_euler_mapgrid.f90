! =========================================================
subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq)
! =========================================================

!     # solve Riemann problems for the 1D Euler equations using Roe's
!     # approximate Riemann solver.

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell
!     # On output, wave contains the waves,
!     #      the speeds,
!     #            amdq the  left-going flux difference  A^- \Delta q
!     #            apdq the right-going flux difference  A^+ \Delta q

!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routine step1, rp is called with ql = qr = q.

    USE Interp_table, ONLY: CO2BLLT_EQUI
    USE var_const, ONLY: e_const,flag_diagno,guessP

    implicit none

    integer :: maxm, meqn, mwaves, mbc, mx, imp, ixy, maux
    double precision ::    ql(1-mbc:maxm+mbc,meqn)
    double precision ::    qr(1-mbc:maxm+mbc,meqn)
    double precision ::     s(1-mbc:maxm+mbc,mwaves)
    double precision ::  wave(1-mbc:maxm+mbc,meqn,mwaves)
    double precision ::  asdq(1-mbc:maxm+mbc,meqn)
    double precision ::  bmasdq(1-mbc:maxm+mbc,meqn)
    double precision ::  bpasdq(1-mbc:maxm+mbc,meqn)
    double precision ::  aux1(1-mbc:maxm+mbc,*)
    double precision ::  aux2(1-mbc:maxm+mbc,*)
    double precision ::  aux3(1-mbc:maxm+mbc,*)

!     # For Roe solver
    double precision :: rhol, ul, vl, el, cl, pl,xl,al,resl,Tl,u2v2l
    double precision :: rhor, ur, vr, er, cr, pr,xr,ar,resr,Tr,u2v2r
    double precision :: u, v, enth, rho, p, e, xm,am,resm,Tm,cm
    double precision :: rhsqrtl, rhsqrtr, rhsq2
    double precision :: uvl, uvr, uv(2),ccm


!     # For mappings
    integer :: meqn2,mwaves2
    parameter(meqn2 = 4, mwaves2 = 3)
    double precision :: ql_state(meqn2), qr_state(meqn2)
    double precision :: q_state(meqn2)
    double precision :: wave_local(mwaves2,meqn2)
    double precision :: s_local(mwaves2), delta(meqn2)
    double precision :: speeds(2,mwaves2)
    double precision :: deltam(meqn2), deltap(meqn2)
    double precision :: area

!     # for mapping
    double precision :: rotm(4), rotp(4), uvm_rot(2),uvp_rot(2)
    integer :: locrot,mcapa,locarea

!     # Miscellaneous
    integer :: i, j, m, mw, i1, ixy1
    logical :: useroe

!     # Problem parameters

    double precision :: dtcom, dxcom, dycom, tcom
    integer :: info

    logical :: in_rpt

    maux  = 7

    useroe = .true.

    call get_aux_locations_t(ixy, mcapa, locrot,locarea)
    do i = 2-mbc,mx+mbc
        i1 = i + imp - 2

        do m = 1,meqn
            ql_state(m) = qr(i-1,m)
            qr_state(m) = ql(i,m)
        enddo

        if (useroe) then
            rhol = ql_state(1)
            rhor = qr_state(1)

            ul = ql_state(2)/rhol
            ur = qr_state(2)/rhor

            vl = ql_state(3)/rhol
            vr = qr_state(3)/rhor

            u2v2l = ul*ul + vl*vl
            u2v2r = ur*ur + vr*vr

            el = ql_state(4)/rhol - e_const - u2v2l
            er = qr_state(4)/rhor - e_const - u2v2r

           CALL CO2BLLT_EQUI(pl, Tl, cl, xl, al, resl,&
     &                         el, 1.0/rhol, 0.d0)
           CALL CO2BLLT_EQUI(pr, Tr, cr, xr, ar, resr,&
     &                         er, 1.0/rhor, 0.d0)

!           # Get Roe averaged values
            rhsqrtl = sqrt(rhol)
            rhsqrtr = sqrt(rhor)
            rhsq2 = rhsqrtl + rhsqrtr
!
            uv(1) = (ul*rhsqrtl + ur*rhsqrtr) / rhsq2
            uv(2) = (vl*rhsqrtl + vr*rhsqrtr) / rhsq2
! ## total enthalpy 
            enth = ( (ql_state(4) + pl)/rhsqrtl + (qr_state(4) + pr)/rhsqrtr )/rhsq2
! ## sound speed average 
            ccm  = (cl*rhsqrtl + cr*rhsqrtr) / rhsq2
        else
        !           # This takes the values needed for the Roe matrix from the
        !           # cell centers in either the  left (imp == 1) or the
        !           # right (imp == 2) cell

            do m = 1,meqn
                if (imp == 1) then
                !                 # Left (minus) cell
                    q_state(m) = ql_state(m)
                else
                !                 # Right (plus) cell
                    q_state(m) = qr_state(m)
                endif
            enddo

            rho   = q_state(1)
            uv(1) = q_state(2)/rho
            uv(2) = q_state(3)/rho
            e     = q_state(4)/rho - 0.5d0*(uv(1)**2 + uv(2)**2)-e_const
            CALL CO2BLLT_EQUI(p, Tm, cm, xm, am, resm,&
     &                         e, 1.0/rho, 0.d0)

            enth = ( q_state(4) + p )/rho
        endif

        do j = 1,2
            uvm_rot(j) = uv(j)
            uvp_rot(j) = uv(j)
            rotm(j) = aux2(i1,locrot+j-1)
            rotp(j) = aux3(i1,locrot+j-1)
        enddo
        call compute_tangent(rotm)
        call compute_tangent(rotp)

        call rotate2(rotm,uv(1),uv(2))
        uvm_rot(1) = uv(1)
        uvm_rot(2) = uv(2)
        call rotate2(rotp,uv(1),uv(2))
        uvp_rot(1) = uv(1)
        uvp_rot(2) = uv(2)

        do m = 1,meqn
            deltap(m) = asdq(i,m)
            deltam(m) = asdq(i,m)
        enddo
        call rotate2(rotm,deltam(2),deltam(3))
        call rotate2(rotp,deltap(2),deltap(3))

    !        # ------------------------------------------
    !        # Solve for minus side
    !        # ------------------------------------------
        ixy1 = 1
        call roe_solver(ixy1,uvm_rot,enth,deltam, &
        wave_local,s_local,info,cm)

        if (info /= 0) then
            write(6,*) 'ixy = ', ixy
            write(6,*) 'imp = ', imp
            write(6,*) 'enth = ', enth
            write(6,*) 'Called from rpt2  (A-DQ)'
            write(6,*) ' '
            write(6,'(24A,24A)') '        Left State      ', &
            '       Right State      '
            write(6,'(24A,24A)') '------------------------', &
            '-----------------------'
            do m = 1,meqn
                write(6,'(2E24.16)') ql_state(m), qr_state(m)
            enddo
            write(6,*) ' '
            stop
        endif

        area = aux2(i1,locarea)
        do mw = 1,mwaves
            call rotate2_tr(rotm,wave_local(mw,2),wave_local(mw,3))
            speeds(1,mw) = area*min(s_local(mw),0.d0)
        enddo

        do m = 1,meqn
            bmasdq(i,m) = 0.d0
            do mw = 1,mwaves
                bmasdq(i,m) = bmasdq(i,m) + speeds(1,mw)*wave_local(mw,m)
            enddo
        enddo

    !        # ------------------------------------------
    !        # Solve for plus side
    !        # ------------------------------------------
        ixy1 = 1
        call roe_solver(ixy1,uvp_rot,enth,deltap, &
        wave_local,s_local,info,cm)

        if (info /= 0) then
            write(6,*) 'ixy = ', ixy
            write(6,*) 'imp = ', imp
            write(6,*) 'enth = ', enth
            write(6,*) 'Called from rpt2  (A+DQ)'
            write(6,*) ' '
            write(6,'(24A,24A)') '        Left State      ', &
            '       Right State      '
            write(6,'(24A,24A)') '------------------------', &
            '-----------------------'
            do m = 1,meqn
                write(6,'(2E24.16)') ql_state(m), qr_state(m)
            enddo
            write(6,*) ' '
            stop
        endif

        area = aux3(i1,locarea)
        do mw = 1,mwaves
            call rotate2_tr(rotp,wave_local(mw,2),wave_local(mw,3))
            speeds(2,mw) = area*max(s_local(mw),0.d0)
        enddo

        do m = 1,meqn
            bpasdq(i,m) = 0.d0
            do mw = 1,mwaves
                bpasdq(i,m) = bpasdq(i,m) + speeds(2,mw)*wave_local(mw,m)
            enddo
        enddo

    enddo !! end of i loop

    return
    end subroutine rpt2
