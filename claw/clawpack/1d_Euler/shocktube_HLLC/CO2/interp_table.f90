MODULE Interp_table
!
!
     USE def_variables     
     USE def_constants
     USE properties
     USE interp_functions
     USE non_linear_solvers
!      
     IMPLICIT NONE
!
     PRIVATE
     PUBLIC  CO2BLLT_EQUI !,CO2BLLT_META, CO2BCLT_EQUI, CO2BCLT_META
!
!
     CONTAINS
!
!===============================================================================
!
SUBROUTINE CO2BLLT_EQUI(p_out, T_out, c_out, x_out, a_out, res, u_in,v_in,p_guess) 
!
!===============================================================================
!
!                       CO2 BILINEAR LOOK-UP TABLES
!
! Input: u_in (the specific internal energy)
!        v_in (the specific volume v_in) 
! Output: p_out (pressure) 
!         T_out (temerature)
!         c_out (speed of sound)
!         x_out (thermodynamic quality for two-phase regions) 
!
!===============================================================================
REAL(pr),INTENT(OUT)  :: p_out, T_out, x_out, a_out, c_out, res
REAL(pr),INTENT(IN)   :: u_in, v_in,p_guess
!
!Local Variables
INTEGER :: i, j, flag_TP, j_sat, i_R, i_L, flag_loca,Niter,exitflag
REAL(pr) :: v_min, v_max, v_sat, v_sat_log, delta,&
&            qual,press,temp,sound,T_guess,out2
REAL(pr) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
REAL(pr) :: vL, uL, vV, uV, pp, ratio, du_dp_x, dv_dp_x
REAL(pr) :: x_check, v_check
!
!=======================================================================================
!
! 1 step: locating the couple (v,e) in the domain physique
!        - HT: if the input value of internal energy is higher than the value of the
!          maximum point on the saturation curve, the point lies in the HT
!          Region.
!
!          If the input value of internal energy is lower than the value which
!          corresponds to the maximum saturation point, it is necessary to 
!          compare the input value of the specific volume with the specific
!          volume of the maximum point of the curve. If the input volume is lower, 
!          the input point is in the left part with respect to the maximum, so the 
!          loop enters the LL_LH section
!
!        - LL_LH:it is necessary to ensure that the point is inside the defined 
!          and tabulated range of validity. The loop stops working if the input 
!          point is outside (internal energy lower than the one referred to the 
!          triple point, pressure higher than the maximum allowed). 
!          If it inside, a comparison between its value of specific internal energy 
!          and the critical one, leads to the choice between LL Region and LH Region.

!        - if the point is instead in the right part with respect to the maximum,
!          a similar reasoning is required. It is checked that the point 
!          belongs to the tabulated range (its pressure cannot be lower than the
!          Triple Point pressure). Then, through a comparison between values of 
!          internal energy and specific volume with respect to the ones at the 
!          saturation, R Region or TP Region are chosen respectively. 
!           
!
!2 step:   choose suitable subroutine to evaluate properties.  
!
!==========================================================================================
!flag_TP = 0
flag_loca = 0
x_out   = 1_pr
a_out   = 1_pr
res     = 0_pr
!
IF (u_in .GT. e_umax) THEN
! In HT
          IF (u_in .GT. u_end) THEN
             STOP '** Out of range. Too high specific internal energy'
          ENDIF
!
! To evaluate the interval on the vertical axis
          delta = y_mesh_HT(2) - y_mesh_HT(1)
          i     = INT((u_in    - y_mesh_HT(1))/delta) + 1 !location indice 
!
          v_min = 0_pr
          v_max = 0_pr
          DO j = 1, ord_spline + 1
             v_min = v_min + spline_left_HT (ord_spline+2-j,i) * u_in**(j-1)
             v_max = v_max + spline_right_HT(ord_spline+2-j,i) * u_in**(j-1)
          ENDDO
!
          IF (v_in .LT. v_min) THEN
             STOP '** Out of range. Pressure higher than 10 MPa in HT'
          ELSEIF (v_in .GT. v_max) THEN
!             STOP '** Out of range. Pressure lower than the triple point pressure in HT'
          ! small pressure region
          flag_loca = 6
          ELSE
          ! Then we are in Region High Temperature (HT)
          flag_loca = 4
          END IF
ELSEIF ((u_in .LE. e_umax) .AND. (v_in .LE. v_umax)) THEN
! In LL or LH
        IF (u_in .LE. e_cr) THEN
! In LL
                IF (u_in .LT. e_tri_L) THEN
                STOP '** Out of range. Too low specific internal energy value in LL'
                ENDIF  
!
! To evaluate the interval on the vertical axis
          delta  = y_mesh_LL(2) - y_mesh_LL(1)
          i      = INT((u_in    - y_mesh_LL(1))/delta) + 1
!
          v_min = 0_pr;          v_sat = 0_pr
          DO j = 1, ord_spline + 1
             v_min = v_min + spline_pmax_LL(ord_spline+2-j,i) * u_in**(j-1)
             v_sat = v_sat + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!             v_sat_log = v_sat_log + spline_Lsat_LL(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat     = 10_pr ** v_sat_log
!print*, spline_Lsat_LL(ord_spline+2-j,i)
          ENDDO
!
!print*,'============================='
!print*, "v_min", v_min
!print*, "v_sat", v_sat
                IF (v_in .LT. v_min) THEN
                STOP '** Out of range. Pressure higher than 10 MPa in LL'
                ELSEIF (v_in .GT. v_sat) THEN
                 ! Region Two-phase (TP)
                        flag_loca = 5  
                ELSE
                 ! Region Left Low (LL) 
                        flag_loca = 1 
                ENDIF
! End in LL and begin in LH
        ELSE
!
!  To evaluate the interval on the vertical axis
          delta = y_mesh_LH(2)- y_mesh_LH(1)
          i     = INT((u_in - y_mesh_LH(1))/delta) + 1
!
          v_min = 0_pr;          v_sat = 0_pr
          v_sat_log = 0_pr
          DO j = 1, ord_spline + 1
             v_min     = v_min     + spline_pmax_LH(ord_spline+2-j,i) *u_in**(j-1)
             v_sat_log = v_sat_log + spline_Lsat_LH(ord_spline+2-j,i) *u_in**(j-1)
             v_sat     = 10_pr ** v_sat_log
          ENDDO
!
          IF (v_in .LT. v_min) THEN
             STOP '** Out of range. Pressure higher than 10 MPa in LH'
          ELSEIF (v_in .GT. v_sat) THEN
          ! Region Two-phase (HT)   
                flag_loca = 5
          ELSE
          ! Region Left High (LH)
                flag_loca = 2
          ENDIF
!End in LH
        ENDIF
! Begin in R or p<5bar
ELSE
! It occurs that:  ((u_in .LE. u_max) .AND. (v_in .GT. v_umax))
!
! To evaluate the interval on the vertical axis
     IF (u_in .LT. e_tri_R) THEN
        x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
        v_check = x_check*v_tri_R + (1_pr-x_check)*v_tri_L
        IF (v_in .GT. v_check) THEN
!       STOP '** Out of range two-phase region'
        ! Region small pressure
        flag_loca = 6
        ENDIF
        ! Region two-phase
        flag_loca = 5
     ELSE      

        delta = y_mesh_R(2) - y_mesh_R(1)
       i     = INT((u_in   - y_mesh_R(1))/delta) + 1
!
       v_sat = 0_pr;       v_max = 0_pr
       DO j = 1, ord_spline + 1
          v_sat = v_sat + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
          v_max = v_max + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
       ENDDO
!
! To evaluate if it is in the single phase domain or in the two-phase
!
!       IF ((u_in .LT. e_umax) .AND. (u_in .GT. e_tri_R) .AND. &
!&           (v_in .LT.v_sat)) THEN
        ! Region two-phase  
!        flag_loca = 5   
!        
!       ENDIF
!
       IF (v_in .GT. v_max) THEN
!          STOP '** Out of range. Pressure lower than the triple point pressure in R'
       ! Small pressure
        flag_loca=6
       ELSEIF (v_in .LT. v_sat) THEN
        ! Region two-phase
          flag_loca = 5
       ELSE
        ! Region Right (R)
          flag_loca = 3
       ENDIF
     ENDIF
!
ENDIF
!print*, "location", flag_loca
!
!
!
SELECT CASE (flag_loca)
!
CASE( 0 )
STOP '** Locating the points in the physical domaion failed '
!
! LL
CASE( 1 )
     CALL Lin_int_Left_Low(T_out, p_out, c_out, u_in, v_in)
     x_out   = 0_pr
     a_out   = 0_pr
! LH
CASE( 2 ) 
     CALL Lin_int_Left_High(T_out, p_out, c_out, u_in, v_in)
! R
CASE( 3 )
     delta        = (x_mesh_max - x_mesh_min)/(MMM_R-1)
     x_mesh_R     = x_mesh_min + (/(i*delta, i=0,MMM_R-1)/)
     CALL Lin_int_Log10(T_out, p_out, c_out, u_in, v_in, NNN_R, MMM_R, x_mesh_R, y_mesh_R, &
&         spline_pmin, spline_Vsat, TTT_R, ppp_R, ccc_R) 
! HT 
CASE( 4 )
     delta       = (x_mesh_max - x_mesh_min)/(MMM_HT-1)
     x_mesh_HT   =  x_mesh_min + (/(i*delta, i=0,MMM_HT-1)/)
     CALL Lin_int_Log10(T_out, p_out, c_out, u_in, v_in, NNN_HT, MMM_HT, x_mesh_HT, y_mesh_HT, &
&                            spline_right_HT, spline_left_HT, TTT_HT, ppp_HT, ccc_HT)
!
! Two-phase domain
CASE( 5 ) 

! initial guess for the saturation pressure (p_guess) 

        CALL New_Rap1D(2, press, qual, res, Niter,&
     &                exitflag, u_in, p_guess, v_in, temp)
!        
        IF (res > 10e-10) THEN
          print*, "res", res, "iter", Niter
          STOP '** Interpolation in two-phase region failed in Interp_table case(5)'
        ENDIF
p_out = press
x_out = qual
T_out = temp
!
! SPEED OF SOUND CALCULATION for the HEM:
!
       delta = saturP(2) - saturP(1)
       j_sat = INT((p_out  - saturP(1))/delta) + 1
!
       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
       DO i = 1, ord_spline
          pp      = p_out**(ord_spline - i)
          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
       ENDDO
!
!       duL_dp = duL_dp * 1e-3_pr              ! (J/kg)/(Pa) 
!       duV_dp = duV_dp * 1e-3_pr              ! (J/kg)/(Pa)
!       dvL_dp = dvL_dp * 1e-6_pr              ! (m3/kg)/(Pa)
!       dvV_dp = dvV_dp * 1e-6_pr              ! (m3/kg)/(Pa)
!
       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr
       DO i = 1, ord_spline+1
          pp = p_out**(i-1)
          vL = vL + vL_psat_spline(ord_spline+2-i, j_sat) * pp
          uL = uL + uL_psat_spline(ord_spline+2-i, j_sat) * pp
          vV = vV + vV_psat_spline(ord_spline+2-i, j_sat) * pp
          uV = uV + uV_psat_spline(ord_spline+2-i, j_sat) * pp
       ENDDO
!
!
       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!
       c_out = SQRT((p_out  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!
       a_out = x_out*vV/v_in
!
       IF (((a_out .GT. 1_pr) .AND. (x_out .GT. 8e-1_pr)) .OR. (x_out .GT. a_out)) THEN  
             a_out = x_out
       ENDIF 
!
 CASE( 6 ) !(perfect gas EoS can be considered and no iteratif process ???)
! Using directly the EoS of span-wagner for small pressure region
        IF (u_in .LT. -123.74e3_pr) THEN
         STOP '** out of range (SOLID) in Interp_table case(6)'
        ENDIF
!
        CALL New_Rap1D(1,temp,out2,res,Niter,&
                        exitflag,u_in,T_guess,v_in,out2)
        IF (res > 10e-10) THEN
          print*, "res", res, "iter", Niter
          STOP '** conversion from internal energy to temperature in small pressure  region failed in Interp_table case(6)'
        ENDIF        
!
        CALL pressure(temp, v_in, press)
        CALL sound_speed(temp, v_in, sound)
!
        T_out = temp
        p_out = press
        c_out = sound
!
        IF ((p_out .GE. 0.5_pr) .OR. (p_out .LT. 0_pr)) THEN
         STOP '** Problem in small pressure region in Interp_table case (6)' 
        ENDIF
!
END SELECT
!
END SUBROUTINE CO2BLLT_EQUI
!
!=============================================================================================
!
!CO2BLLT_META, CO2BCLT_EQUI, CO2BCLT_META
!
!============================================================================================
!
!
!
!
END MODULE Interp_table 
