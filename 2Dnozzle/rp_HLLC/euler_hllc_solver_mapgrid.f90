    subroutine hllc_solver(pl,prr,ul,ur,vl,vr,cl,cr,rhol,rhor,&
&                          ql_state,qr_state,wave_local,s_local,info)
    USE def_constants, ONLY: pr
    implicit none

    REAL(pr) :: pl,prr,ul,ur,vl,vr,cl,cr,rhol,rhor
    REAL(pr) :: wave_local(3,4), s_local(3)
    REAL(pr) :: sl,sr,sm
    REAL(pr) :: qml(4),qmr(4)
    REAL(pr) :: ql_state(4),qr_state(4)
    integer :: m, p, mu, mv, ixy, i, j, k, info

    info = 0
!    print*, 'HLLC solver called'
!     # These are used even in the mapped case, but ixy is
!     # always set to 1;  rotation matrix rot does switches
!     # coordinates for us.
      sl     = min(ul-cl, ur-cr)
      sr     = max(ul+cl, ur+cr)
      sm     = (prr-pl + rhol*ul*(sl- ul) - rhor*ur*(sr- ur))&
&                  / (rhol*   (sl- ul) - rhor*   (sr- ur))
!      print*, sl,sr,sm
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
      qmr(4) = (qr_state(4)/rhor + (sm - ur)*(sm + prr/     &
&                        (rhor*(sr - ur))))
!         qmr(1:4) = qmr(1:4) * rhor*(sr - ur)/(sr - sm)
!  # Compute the waves.
!
!      print*, qml(1),qml(2),qml(3),qml(4)
      !print*, qml(1),qml(2),qml(3),qml(4)
      DO m =1,4
         qmr(m) = qmr(m) * rhor*(sr - ur)/(sr - sm)
         qml(m) = qml(m) * rhol*(sl - ul)/(sl - sm)
      ENDDO
!     # 1-wave (rarefaction or shock)
    wave_local(1,1)  = qml(1) - ql_state(1)
    wave_local(1,2) = qml(2) - ql_state(2)
    wave_local(1,3) = qml(3) - ql_state(3)
    wave_local(1,4)  = qml(4) - ql_state(4)
    s_local(1) = sl

!     # contact discontinuity
    wave_local(2,1)  = qmr(1)  - qml(1)
    wave_local(2,2) = qmr(2)  - qml(2)
    wave_local(2,3) = qmr(3)  - qml(3)
    wave_local(2,4)  = qmr(4)  - qml(4)
    s_local(2) = sm

!     # 3-wave (rarefaction or shock)
    wave_local(3,1)  = qr_state(1) - qmr(1)
    wave_local(3,2) = qr_state(2) - qmr(2)
    wave_local(3,3) = qr_state(3) - qmr(3)
    wave_local(3,4)  = qr_state(4) - qmr(4)
    s_local(3) = sr

    end subroutine hllc_solver
