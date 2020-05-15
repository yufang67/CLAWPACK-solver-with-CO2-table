! ===================================================================
!
!                      Saturation_Curve 
!
! ===================================================================
! This subroutine create a grid in the 2-phase region by discretizing P_sat
! from P_tri to P_crit to obtain a table of the saturation properties:
!
!                 P_sat, V_v, V_l, e_v, e_l, T_sat 
!
! The spline coefficients are also computed, so for given values of 
! saturation pressure, coefficients of spline are ble to express 
!                   u(p_sat), v(p_sat) , T(p_sat)
!
! ===================================================================
SUBROUTINE saturation_curve()
       
     USE def_constants
     USE def_variables
     USE non_linear_solvers
     USE properties
     USE grid_functions 
        IMPLICIT NONE        
!
        INTEGER  :: i, j, exitflag,Niter
!
        REAL(pr)  :: delta, delta_sat,v_v,v_l,Tsat,u_l,u_v,&
 &                   resnorm, guess1,guess2,guess3,in_2
!
        REAL(pr), DIMENSION (NNN_sat_TP) :: vL, vV, uL, uV,saturT
!
!
!
        delta      = (P_cr - P_tri) / (NNN_TP-1)
        saturP     = P_tri + (/(i*delta, i=0, NNN_TP-1)/)
        delta_sat  = (P_cr - P_tri) / (NNN_sat_TP-1)
        saturP_sat = P_tri + (/(i*delta_sat, i=0,NNN_sat_TP-1)/)
!
! points on saturation curve 
!
        guess1 = 1_pr/rho_tri_L
        guess2 = 1_pr/rho_tri_R
        guess3 = T_tri
!
        DO i = 1, NNN_sat_TP-1

         CALL New_Rap3D(3,v_l , v_v, Tsat, &
     &   resnorm, Niter, exitflag, saturP_sat(i),in_2, guess1, guess2,guess3)
        
        IF (resnorm > 1e-9_pr) THEN
        print*, "saturation curve", resnorm, i
        STOP
        ENDIF
!
!        print*, v_l,v_v
!        print*, u_l,u_v
!        print*," "
         CALL inter_energy(Tsat,v_v,u_v)
         CALL inter_energy(Tsat,v_l,u_l)
           vL(i)     = v_l
           vV(i)     = v_v
           uL(i)     = u_l
           uV(i)     = u_v
           saturT(i) = Tsat
!
           guess1 = v_l
           guess2 = v_v
           guess3 = Tsat
        ENDDO
!
        vL(NNN_sat_TP)     = 1_pr/rho_cr
        vV(NNN_sat_TP)     = 1_pr/rho_cr
        uL(NNN_sat_TP)     = e_cr
        uV(NNN_sat_TP)     = e_cr
        saturT(NNN_sat_TP) = T_cr
! intervals on saturation curve
!
        DO i = 1, NNN_TP-1
           j = 3 * (i-1) + 1
             CALL polyfit(vL_psat_spline(:,i), saturP_sat (j:j+3), vL(j:j+3),ord_spline)
             CALL polyfit(vV_psat_spline(:,i), saturP_sat (j:j+3), vV(j:j+3),ord_spline)
             CALL polyfit(uL_psat_spline(:,i), saturP_sat (j:j+3), uL(j:j+3),ord_spline)
             CALL polyfit(Tsat_psat_spline(:,i), saturP_sat(j:j+3),saturT(j:j+3), ord_spline)
             CALL polyfit(uV_psat_spline(:,i), saturP_sat (j:j+3), uV(j:j+3),ord_spline)
        ENDDO
print*,'-------------------------------------------------------------------------'
print*,'                  TP SATURATION CURVE FINISH                             '
print*,'-------------------------------------------------------------------------'
!
!
      END SUBROUTINE saturation_curve
