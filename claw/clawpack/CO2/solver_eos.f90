MODULE solver_eos
!
!
!USE properties
!USE peng_robinson
!USE interp_table
USE var_const, ONLY:  e_const
USE def_constants, ONLY: pr,ord_spline
USE def_variables, ONLY: vL_psat_spline, vV_psat_spline,&
                        &uL_psat_spline, uV_psat_spline,&
                        &Tsat_psat_spline, saturP

!
IMPLICIT NONE
PUBLIC  eos_1d
CONTAINS


!===============================================================================
          SUBROUTINE eos_1d(MODE, out_1, out_2, resnorm, Niter,&
     &                             exitflag, GGG, X0, in_1, out3)
!===============================================================================
!         1D        inverse the internal energy to find the temperature
!------------------------------------------------------------------------------
!
!  Input :
!  -------
!     MODE        =  2 Peng-robinson, 4 span-wagner, 5 (P,T)--> v and e
!     GGG         =  e
!     in_1        =  v
!     X0          =  initial guess of T
!
!
!  Output :
!  -------
!
!     out_1       =  final T
!     out_2       =  idem
!     out3       =  idem
!
!-------------------------------------------------------------------
!                                         M. De Lorenzo      03/2016
!-------------------------------------------------------------------
!
!
      IMPLICIT NONE
!
      INTEGER             ::  MAX_ITER
      INTEGER, INTENT(IN) ::  MODE
      INTEGER, INTENT(OUT) :: exitflag, Niter
      REAL(pr), INTENT(IN)   :: GGG, X0, in_1
      REAL(pr), INTENT(OUT)  :: out_1, out_2, resnorm, out3
      REAL(pr)  ::  TOL_X, TOL_FUN, ALPHA, MIN_LAMBDA, MAX_LAMBDA, Jstar,&
     &  lambda, slope, fffold, lambda_min, fff, fff2, A_l, aa, bb, Dx_J,&
     &  discriminant, lambda2, lambda_OLD, XXX, F, delta, XXX_forw, J0,&
     &  XXX_back, F_forw, F_back, dF, TYPX, dx_star, dx, XXX_old, g, J
      REAL(pr), DIMENSION (2,2)  ::  B_l
      REAL(pr), DIMENSION (2,1)  ::  C_l, aabb
!
!
      XXX = X0                        ! Initial value
!
!   Options
!
      TOL_X      = 1d-20            ! relative max step tolerance
      TOL_FUN    = 1d-12            ! function tolerance
      MAX_ITER   = 500              ! max number of iterations
      ALPHA      = 1d-4             ! criteria for decrease
      MIN_LAMBDA = 1d-1             ! min lambda
      MAX_LAMBDA = 5d-1
      TYPX       = abs(XXX)           ! x scaling value, remove zeros
!
! ---- Function F definition
!      F       = GGG_input - GGG
!
      call New_Rap1D(MODE, F, out_2, XXX, GGG, in_1, out3)
!
!
! --------- Jacobian estimation-----------------
!
      Dx_J       = 1d-6         ! finite difference delta
      delta      = Dx_J
      XXX_forw   = XXX + delta
      XXX_back   = XXX - delta
      call New_Rap1D(MODE, F_forw, out_2, XXX_forw, GGG, in_1, out3)
      call New_Rap1D(MODE, F_back, out_2, XXX_back, GGG, in_1, out3)
      dF       = F_forw - F_back    ! delta F
      J        = dF/Dx_J/2d0        ! derivatives dF/d_n
      J0       = 1d0 / TYPX         ! Jacobian scaling factor
      Jstar    = J/J0;              ! scale Jacobian
      resnorm  = abs(F)             ! Norm of residues
      exitflag = 1
!
!
!--- SOLVER -------------------
!
      Niter  = 0
      lambda = 1d0                ! backtracking
      lambda_OLD = lambda
!
!
!
1       DO WHILE ( ((resnorm>TOL_FUN) .OR. (lambda<1d0)) .AND. &
     &             (exitflag>=0) .AND. (Niter<=MAX_ITER))
             Niter = Niter + 1
!
!--------- Newton-Raphson solver
!
             IF (lambda .eq. 1.0d0) THEN
               dx_star = -F/Jstar        ! calculate Newton step
               dx      = dx_star*TYPX
               g       = F*Jstar
               slope   = g*dx_star
               fffold  = F*F
               XXX_old = XXX
               lambda_min = TOL_X/(ABS(dx)/ABS(XXX_old))
            ENDIF
!
!--------- Check about proximity of XXX and XXX_OLD
!
             IF (lambda < lambda_min) THEN
                exitflag = 2;
                !print*,'XXX is too close to XXX_OLD'
                !EXIT ! OUT NewRap
             ENDIF
!
!--------- Eventually, backtracking of New-Rap step
!
             XXX = XXX_old + dx*lambda
             call New_Rap1D(MODE, F, out_2, XXX, GGG, in_1, out3)
!
!--------- Jacobian estimation
!
             XXX_forw   = XXX + delta
             XXX_back   = XXX - delta
             call New_Rap1D(MODE, F_forw, out_2, XXX_forw, GGG, in_1,&
     & out3)
             call New_Rap1D(MODE, F_back, out_2, XXX_back, GGG, in_1,&
     & out3)
             dF     = F_forw - F_back    ! delta F
             J      = dF/Dx_J/2d0        ! derivatives dF/d_n
             Jstar  = J/J0               ! scale Jacobian
             fff    = F*F                ! f to be minimized
!
!--------- Optimization technique (minimizing fff = SUM(F*F) )
!
             IF (fff > fffold + ALPHA*lambda*slope) THEN
                 IF (lambda .eq. 1d0) THEN
                     lambda = -slope/2d0/(fff-fffold-slope) ! calculate lambda
                 ELSE
                     A_l      =  1d0/(lambda_OLD - lambda2)
                     B_l(1,1) =  1d0/lambda_OLD/lambda_OLD
                     B_l(1,2) = -1d0/lambda2/lambda2
                     B_l(2,1) = -lambda2/lambda_OLD/lambda_OLD
                     B_l(2,2) =  lambda_OLD/lambda2/lambda2
                     C_l(1,1) =  fff-fffold-lambda_OLD*slope
                     C_l(2,1) =  fff2-fffold-lambda2*slope
                     !
                     aabb = A_l* MATMUL(B_l,C_l)
                     aa   = aabb(1,1)
                     bb   = aabb(2,1)
                     IF (aa .EQ. 0.0d0) THEN
                         lambda = -slope/2d0/bb;
                     ELSE
                         discriminant = bb**2d0 - 3d0*aa*slope;
                         IF (discriminant < 0d0) THEN
                             lambda = MAX_LAMBDA*lambda_OLD;
                         ELSEIF (bb <= 0d0) THEN
                             lambda = (-bb + sqrt(discriminant))/3d0/aa
                         ELSE
                             lambda = -slope/(bb + sqrt(discriminant))
                         ENDIF
                     ENDIF
                 lambda = min(lambda,MAX_LAMBDA*lambda_OLD)     ! minimum step length
                 ENDIF
             ELSE
                 lambda = 1d0; ! fraction of Newton step
             ENDIF
             IF (lambda < 1d0) THEN
                 lambda2 = lambda_OLD;
                 fff2    = fff;       ! save 2nd most previous value
                 lambda  = max(lambda,MIN_LAMBDA*lambda_OLD) ! minimum step length
                 GO TO 1
             ENDIF
!
             resnorm = abs(F)
!
      END DO
!
      out_1 = XXX
!

      END SUBROUTINE eos_1d
!
!
!===============================================================================
      SUBROUTINE New_Rap1D(MODE, F, out_2, XXX, GGG, in_1, out_3)
!===============================================================================
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: MODE
          REAL(pr), INTENT(IN)  :: in_1, XXX, GGG
          REAL(pr), INTENT(OUT) :: F, out_2,out_3
!
!LOCAL
          INTEGER :: i,j
          REAL(pr) :: temp, e, v,sound,a_out,res,p_guess,eint,vin
          REAL(pr) :: delta_p,press,e_l,e_v,V_v,V_l,qual
          
         

! (v,e) --> T: v=in_1, e=GGG, T=XXX
!          v = in_1
!          temp = XXX
         
          IF (MODE .EQ. 4) THEN
!                v = in_1
!                temp = XXX
!                       
!                CALL inter_energy(temp, v, e)
!                F = (e - GGG) * 1d-3
!                out_2 = 0.0d0
!                out_3 = 0.0d0
           ELSEIF (MODE .EQ. 1) THEN
! In two-phase region, (e,v) --> (psat,x,Tsat) 
!       
        press = XXX
        vin = in_1
        eint = GGG + 500e3_pr
        delta_p = saturP(2) - saturP(1)
        i = INT((press - saturP(1))/delta_p) + 1
!
        
        e_l = 0_pr
        e_v = 0_pr
        v_l = 0_pr
        v_v = 0_pr
        temp = 0_pr
!
        DO j = 1, ord_spline + 1
        e_l = e_l + uL_psat_spline(ord_spline+2-j,i) * press**(j-1)
        e_v = e_v + uV_psat_spline(ord_spline+2-j,i) * press**(j-1)
        v_l = v_l + vL_psat_spline(ord_spline+2-j,i) * press**(j-1)
        v_v = v_v + vV_psat_spline(ord_spline+2-j,i) * press**(j-1)
        temp = temp + Tsat_psat_spline(ord_spline+2-j,i) * press**(j-1)
!
        ENDDO
!        print*, 'iter process',v,v_l, v_v
        out_3 = temp
        qual = (vin - v_l)/(v_v-v_l)
!
        IF( (qual /= qual) .OR. (qual < 0.0)) then 
          print*, 'eos_1d qual problem', vin, v_l,v_v,press
        ENDIF
!
        out_2 = qual
        e_v = e_v + 500e3_pr
        e_l = e_l + 500e3_pr
        F = (eint - qual*e_v - (1_pr - qual)*e_l) * 1e-3_pr
        IF (F /= F) then
          print*, 'eos_1d F problem', eint,e_v,e_l
        ENDIF
!        print*, 'qual', qual,'e',e,qual*e_v + (1_pr - qual)*e_l
!               
          ELSEIF (MODE .EQ. 2) THEN
!               
!                v = in_1
!                temp = XXX              
 
!                CALL interenergy_pr(temp, v, e)
!                F = (e - GGG) * 1d-3
!                out_2 = 0.0d0
!                out_3 = 0.0d0
!
          ELSEIF (MODE .EQ. 5) THEN
! (T,p) --> rho,e, temp=in_1,
 !               temp = in_1     
 !               v = XXX
 !               CALL pressure(temp,v,press)
 !               F = (press - GGG) * 10d-6 
 !              
 !               CALL inter_energy(temp,v,e)  
 !               out_2 = e
          ELSEIF (MODE .EQ. 3) THEN
! (v,p) --> e  for outletBC using LOOKUP-table
!                v = in_1
!                e = XXX
!                print*, 1.0/v, e
!                       
!                CALL CO2BLLT_EQUI(press,temp,sound,&
!      &              qual, a_out, res, e,v,5d6)
!                F = (press - GGG) * 1d-6
!                guessp_BC(i_bc) = press
!
!                out_2 = 0.0d0
!                out_3 = 0.0d0                
!
          ELSEIF (MODE .EQ. 6) THEN
!(e,p)--> v  for outlet BC
!                e = in_1
!                v = XXX
!                       
!                guessp_BC(i_bc) = (guessp_BC(i_bc)*2.0+guessP(100,i_bc))/3.0
!                print*,'guessp_BC', guessp_BC(i_bc)
!                p_guess = (3.0*guessp_1d(i_bc) + guessp_1d(i_bc-1))/4.0
!                CALL CO2BLLT_EQUI(press,temp,sound,&
!      &              qual, a_out, res, e,v,p_guess)
!                F = (press - GGG) * 1d-6
!                print*, 'press_outlet',press
!                guessp_BC(i_bc) = press
!                guessp_1d(i_bc) = press
!
!                out_2 = 0.0d0
!                out_3 = 0.0d0
          ELSE

                print*,'MODE of EOS unknown'
                STOP
          ENDIF          
!
      END SUBROUTINE New_Rap1D
END MODULE solver_eos
