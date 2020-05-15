MODULE Interp_table
!
!
     USE def_variables     
     USE def_constants
     USE properties
     USE interp_functions
     USE non_linear_solvers
     USE var_const, ONLY: flag_diagno,flag_out
     USE solver_eos, ONLY:eos_1d
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
! NOTE: TP region tabulated and mix with NewtonRaph to optimise the errors except 
!       a small region close to the critical point 27/11/2017.
!       The boundary between LH and HT is optimised. 25/11/2017     
!===============================================================================
REAL(pr),INTENT(OUT)  :: p_out, T_out, x_out, a_out, c_out, res
REAL(pr),INTENT(IN)   :: p_guess
!
!Local Variables
INTEGER :: i, j, flag_TP, j_sat, i_R, i_L, flag_loca,Niter,exitflag
REAL(pr) :: v_min, v_max, v_sat, v_sat_log, delta,&
&            qual,press,temp,sound,T_guess,out2
REAL(pr) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
REAL(pr) :: vL, uL, vV, uV, pp, ratio, du_dp_x, dv_dp_x
REAL(pr) :: x_check, v_check
REAL(pr) :: gamma_pg, e_pg, v_in, u_in
REAL(pr) :: resnorm,out3,T_gs, p_gs, c_gs,dummy
REAL(pr) :: eint, vin, pguess,res2
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
flag_TP = 1
flag_loca = 0
x_out   = 1_pr
a_out   = 1_pr
res     = 0_pr
!
IF (u_in .GT. e_umax) THEN
! 
          IF (u_in .GT. u_end) THEN
             print*, 'diagno=', flag_diagno
             STOP '** Out of range. Too high specific internal energy'
          ENDIF
!
! To evaluate the interval on the vertical axis, BC of HT
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
             print*, 'diagno=', flag_diagno
             STOP '** Out of range. Pressure higher than 50 MPa in HT'
          ELSEIF (v_in .GT. v_max) THEN
! In LP
             flag_loca = 6
          ELSE
! In HT
             flag_loca = 4
          END IF
ELSEIF ((u_in .LE. e_umax) .AND. (v_in .LE. v_umax)) THEN 
! Left part LL or LH or TP
       IF (u_in .LE. e_cr) THEN
! In LL
           IF (u_in .LT. e_tri_L) THEN
              print*, 'diagno=', flag_diagno
              STOP '** Out of range. Too low specific internal energy value in LL'
           ENDIF  
!
! To evaluate the interval on the vertical axis, BC of LL
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
          IF (v_in .LT. v_min) THEN
             print*, 'diagno=', flag_diagno
             STOP '** Out of range. Pressure higher than 50 MPa in LL'
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
! To evaluate the interval on the vertical axis, BC of LH
          delta = y_mesh_LH(2)- y_mesh_LH(1)
          i     = INT((u_in - y_mesh_LH(1))/delta) + 1
!
          v_min = 0_pr;          v_sat = 0_pr
          v_sat_log = 0_pr
          DO j = 1, ord_spline + 1
             v_min     = v_min     + spline_pmax_LH(ord_spline+2-j,i) *u_in**(j-1)
             v_sat     = v_sat     + spline_Lsat_LH(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat_log = v_sat_log + spline_Lsat_LH(ord_spline+2-j,i) *u_in**(j-1)
!             v_sat     = 10_pr ** v_sat_log
          ENDDO
!
          IF (v_in .LT. v_min) THEN
             print*, 'diagno=', flag_diagno
             STOP '** Out of range. Pressure higher than 10 MPa in LH'
          ELSEIF (v_in .GT. v_sat) THEN
! Region Two-phase (HT)   
                flag_loca = 5
          ELSE
! Region Left High (LH)
                flag_loca = 2
          ENDIF
       ENDIF
! Right part R or TP or LP or solid(not considered yet) 
ELSE
! 
       IF (u_in .LT. e_tri_R) THEN
          x_check = (u_in - e_tri_L)/(e_tri_R - e_tri_L)
          v_check = x_check*v_tri_R + (1_pr-x_check)*v_tri_L
          IF (v_in .GT. v_check) THEN
! In Solid, we put it back to TP
             v_in = v_check - 1e-5_pr
             flag_loca = 5
!             f_out = 1_pr
          ENDIF
! In TP
          flag_loca = 5
       ELSE      
! BC of R
          delta = y_mesh_R(2) - y_mesh_R(1)
          i     = INT((u_in   - y_mesh_R(1))/delta) + 1
!
          v_sat = 0_pr;       v_max = 0_pr
          DO j = 1, ord_spline + 1
            v_sat = v_sat + spline_Vsat(ord_spline+2-j,i) * u_in**(j-1)
            v_max = v_max + spline_pmin(ord_spline+2-j,i) * u_in**(j-1)
          ENDDO
!
          IF (v_in .GT. v_max) THEN
! In LP
            flag_loca=6
          ELSEIF (v_in .LT. v_sat) THEN
! In TP
            flag_loca = 5
          ELSE
! In Region Right (R)
            flag_loca = 3
          ENDIF
       ENDIF
!
ENDIF
!print*, "location", flag_loca
!
flag_out = flag_loca
!
SELECT CASE (flag_loca)
!
CASE( 0 )
STOP '** Locating the points in the physical domaion failed '
!
!## LL ##
CASE( 1 )
     CALL Lin_int_Left_Low(T_out, p_out, c_out, u_in, v_in)
     IF ( (u_in >-200.0e3_pr) .AND. (v_in >1.9e-3_pr) ) THEN
        CALL sound_speed(T_out,v_in,c_out)
        
     ENDIF 
     x_out   = 0_pr
     a_out   = 0_pr
!## LH ##
CASE( 2 )
!#######################################
   IF ( (u_in > e_F) .OR. ((u_in < e_G) .AND. (v_in>v_G)) ) THEN
    CALL Lin_int_Left_High(T_gs, p_gs, c_gs, u_in, v_in)
   
    CALL  New_Rap1D(1, T_out, out2, resnorm, Niter,&
     &              exitflag, u_in,T_gs,v_in, out3) 
    CALL  pressure(T_out,v_in,p_out)
    CALL  sound_speed(T_out,v_in,c_out)
!
   ELSE 
     CALL Lin_int_Left_High(T_out, p_out, c_out, u_in, v_in)
   ENDIF
!########################################
!## R ##
CASE( 3 )
     delta        = (x_mesh_max - x_mesh_min)/(MMM_R-1)
     x_mesh_R     = x_mesh_min + (/(i*delta, i=0,MMM_R-1)/)
     CALL Lin_int(T_out, p_out, c_out, u_in, v_in, NNN_R, MMM_R, x_mesh_R, y_mesh_R, &
&         spline_pmin, spline_Vsat, TTT_R, ppp_R, ccc_R) 
!## HT ## 
CASE( 4 )
     delta       = (x_mesh_max - x_mesh_min)/(MMM_HT-1)
     x_mesh_HT   =  x_mesh_min + (/(i*delta, i=0,MMM_HT-1)/)
     CALL Lin_int_Log10(T_out, p_out, c_out, u_in, v_in, NNN_HT, MMM_HT, x_mesh_HT, y_mesh_HT, &
&         spline_right_HT, spline_left_HT, TTT_HT, ppp_HT, ccc_HT)
!## TP ##
CASE( 5 )
!----------------- Spline Tsat and Psat construction flag_TP=0-------------------------
  IF (flag_TP==0) THEN
! initial guess for the saturation pressure (p_guess) 

!        CALL New_Rap1D(2, press, qual, res, Niter,&
!     &                exitflag, u_in, p_guess, v_in, temp)
!        
!        IF (res > 1e-5) THEN
!          print*, "res IN Interp_table", res, "u,v",u_in,v_in,&
!     &            'T=',temp,'p=',press,'x=',qual
!          STOP '** Interpolation in two-phase region failed in Interp_table case(5)'
!        ENDIF
!p_out = press
!x_out = qual
!T_out = temp
!
! SPEED OF SOUND CALCULATION for the HEM:
!
!       delta = saturP(2) - saturP(1)
!       j_sat = INT((p_out  - saturP(1))/delta) + 1
!
!       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
!       DO i = 1, ord_spline
!          pp      = p_out**(ord_spline - i)
!          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
!          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
!          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
!          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
!       ENDDO
!
!
!       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr
!       DO i = 1, ord_spline+1
!         pp = p_out**(i-1)
!          vL = vL + vL_psat_spline(ord_spline+2-i, j_sat) * pp
!          uL = uL + uL_psat_spline(ord_spline+2-i, j_sat) * pp
!          vV = vV + vV_psat_spline(ord_spline+2-i, j_sat) * pp
!          uV = uV + uV_psat_spline(ord_spline+2-i, j_sat) * pp
!       ENDDO
!
!       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!!
!       c_out = SQRT((p_out  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!!
!       a_out = x_out*vV/v_in
!!
!       IF (((a_out .GT. 1_pr) .AND. (x_out .GT. 8e-1_pr)) .OR. (x_out .GT. a_out)) THEN  
!             a_out = x_out
!       ENDIF 
!------------------ Look-up table in TPL, TPM and TPH, flag_TP=1-----------------------
  ELSEIF (flag_TP==1) THEN
!#########TPH
   IF (u_in .GE. e_tri_R) THEN  
! 
         flag_out = 6
!
         delta       = (x_mesh_max - x_mesh_min)/(MMM_TPH-1)
         x_mesh_TPH  =  x_mesh_min + (/(i*delta, i=0,MMM_TPH-1)/)
         CALL Lin_TP(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPH, MMM_TPH,&
&             x_mesh_TPH, y_mesh_TPH,spline_Vsat, spline_left_TPH,&
&             TTT_TPH, ppp_TPH, ccc_TPH, xxx_TPH)
!
         CALL New_Rap1D(2, press, qual, res, Niter,&
&                       exitflag, u_in, p_out, v_in, temp)
!
         p_out = press
!         T_out = temp
         x_out = qual
!
         delta = saturP(2) - saturP(1)
         j_sat = INT((p_out  - saturP(1))/delta) + 1

       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
       DO i = 1, ord_spline
          pp      = p_out**(ord_spline - i)
          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
       ENDDO
!
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
!##########TPL
   ELSEIF (u_in .LE. e_cr) THEN
!
         flag_out = 8
!         print*, 'in table in', u_in, v_in
         delta       = (x_mesh_max - x_mesh_min)/(MMM_TPL-1)
         x_mesh_TPL  =  x_mesh_min + (/(i*delta, i=0,MMM_TPL-1)/)
!         CALL Lin_TP(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPL, MMM_TPL,&
!&             x_mesh_TPL,y_mesh_TPL,spline_right_TPL, spline_Lsat_LL,&
!&             TTT_TPL, ppp_TPL, ccc_TPL, xxx_TPL)
	 CALL Lin_TP_Log10(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPL, MMM_TPL,&
&             x_mesh_TPL,y_mesh_TPL,spline_right_TPL, spline_Lsat_LL,&
&             TTT_TPL, ppp_TPL, ccc_TPL, xxx_TPL)
	 CALL satprop (3, p_out, dummy, vV, vL, uV, uL)

!        
!         print*, 'in table out', T_out, p_out,c_out,'in',u_in,v_in
!      IF( (u_in < -210e3_pr) .OR.(v_in > 2.34e-3_pr) )THEN 
!         print*, 'in table iter start', p_out, u_in, v_in
!           pguess = p_out 
!           eint   = u_in
!           vin    = v_in
!         CALL eos_1d(1, press, qual, res2, Niter, exitflag,eint, pguess, vin, temp)
!         CALL New_Rap1D(2, press, qual, res, Niter,&
!&                          exitflag, u_in, p_out, v_in, temp)
!
!             print*, 'in table iter out', press, qual
!         IF (res2 > 1_pr) THEN
!            print*, 'TPL interpolation guess prob', res2 
!            print*, 'u,v',u_in,v_in, 'p_iter',press, 'p_tb', p_out
!         ENDIF
!         p_out = press
!         x_out = qual
!
!         delta = saturP(2) - saturP(1)
!         j_sat = INT((p_out  - saturP(1))/delta) + 1
!
!       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
!       DO i = 1, ord_spline
!          pp      = p_out**(ord_spline - i)
!          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
!          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
!          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
!          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
!       ENDDO
!
!
!       vL = 0_pr; uL = 0_pr;  vV = 0_pr;  uV = 0_pr
!       DO i = 1, ord_spline+1
!          pp = p_out**(i-1)
!          vL = vL + vL_psat_spline(ord_spline+2-i, j_sat) * pp
!          uL = uL + uL_psat_spline(ord_spline+2-i, j_sat) * pp
!          vV = vV + vV_psat_spline(ord_spline+2-i, j_sat) * pp
!          uV = uV + uV_psat_spline(ord_spline+2-i, j_sat) * pp
!       ENDDO
!
       x_out   = (v_in - vL) / (vV - vL)
!       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!
!       c_out = SQRT((p_out  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!         
       a_out = x_out*vV/v_in
!
         IF (((a_out .GT. 1_pr) .AND. (x_out .GT. 8e-1_pr)) .OR. (x_out .GT. a_out)) THEN
            a_out = x_out
         ENDIF
!      ENDIF
!###############TPM
   ELSE  
!
         flag_out = 7
         delta       = (x_mesh_max - x_mesh_min)/(MMM_TPM-1)
         x_mesh_TPM  =  x_mesh_min + (/(i*delta, i=0,MMM_TPM-1)/)
         CALL Lin_TP_Log10(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPM, MMM_TPM,&
&             x_mesh_TPM,y_mesh_TPM,spline_right_TPM, spline_left_TPM,&
&             TTT_TPM, ppp_TPM, ccc_TPM, xxx_TPM)
!
!print*, 'p',p_out
       IF ( (u_in > -183e3_pr) .OR. (v_in >2.3408e-3_pr) ) THEN    
          pguess = p_out
          eint   = u_in
          vin    = v_in
          CALL eos_1d(1, press, qual, res2, Niter, exitflag,eint, pguess, vin, temp)
!print*, 'after eos_1d', press
!         CALL New_Rap1D(2, press, qual, res, Niter,&
!     &                exitflag, u_in, p_out, v_in, temp)
!
         IF (res > 1_pr) THEN
            print*, 'TPM interpolation guess prob', res
            print*, 'u,v',u_in,v_in, 'p_iter',press, 'p_tb', p_out
         ENDIF
         p_out = press
         x_out = qual
         IF (x_out < 1e-8_pr) THEN
             x_out=0.0
         ENDIF
!
         delta = saturP(2) - saturP(1)
         j_sat = INT((p_out  - saturP(1))/delta) + 1

       duL_dp = 0_pr;  duV_dp = 0_pr;  dvL_dp = 0_pr;  dvV_dp = 0_pr
       DO i = 1, ord_spline
          pp      = p_out**(ord_spline - i)
          duL_dp  = duL_dp + (ord_spline+1-i) * uL_psat_spline(i,j_sat) *pp
          duV_dp  = duV_dp + (ord_spline+1-i) * uV_psat_spline(i,j_sat) *pp
          dvL_dp  = dvL_dp + (ord_spline+1-i) * vL_psat_spline(i,j_sat) *pp
          dvV_dp  = dvV_dp + (ord_spline+1-i) * vV_psat_spline(i,j_sat) *pp
       ENDDO
!
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
      ENDIF    
   ENDIF
  ENDIF
!## LP ##  perfect gas OR  span-wagner (iterative)
CASE( 6 )
!        
        IF (u_in .LT. -123.74e3_pr) THEN
! In Solid phase, but we push it back to LP
         u_in  = -123.74e3_pr + 1e3_pr
!         f_out = 2_pr
         Print*,'** out of range (SOLID) in Interp_table case(6)'
        ENDIF
! Using directly the EoS of span-wagner for small pressure region
!
!        CALL New_Rap1D(1,temp,out2,res,Niter,&
!                        exitflag,u_in,T_guess,v_in,out2)
!        IF (res > 10e-10) THEN
!          print*, "res IN Interp_table", res, "iter", Niter
!          STOP '** conversion from internal energy to temperature in small pressure  region failed in Interp_table case(6)'
!        ENDIF        
!
!        CALL pressure(temp, v_in, press)
!        CALL sound_speed(temp, v_in, sound)
!
!        T_out = temp
!        p_out = press
!        c_out = sound
!
!        IF ((p_out .GE. 0.5_pr) .OR. (p_out .LT. 0_pr)) THEN
!         STOP '** Problem in small pressure region in Interp_table case (6)' 
!        ENDIF
!--------------  perfect gas formulation--------------------------------------------- 
         gamma_pg = 1.313_pr
         e_pg     = 236.0294e3_pr
         p_out    = (gamma_pg-1_pr)*(u_in+e_pg)/v_in
         c_out    = sqrt(gamma_pg*p_out*v_in)
         T_out    = p_out*v_in/R
         x_out    = 1_pr
         a_out    = 1_pr
         IF ( (p_out .LT. 1e-2_pr) .OR. (c_out .LT. 1e-2_pr ) .OR. (T_out .LT. 1e-2_pr)) THEN
            Print*, 'p=',p_out, 'T=', T_out, 'c', c_out
            STOP '** Problem Perfect gas formulation in Interp_table case (6)' 
         ENDIF
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
END MODULE Interp_table 
