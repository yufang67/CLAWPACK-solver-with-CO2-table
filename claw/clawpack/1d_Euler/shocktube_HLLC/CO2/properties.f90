!---------------------------------------------------------------------------------------------------
! All the formulations  have been taken from the original
! article,the Span-Wagner EoS for CO2.
! The EoS has been published in J. Phys. Chem. Ref. Data, Vol.25,
! pp. 1509-1596, 1996.
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! MODULE   properties
! @brief   Compute thermodynamic properties  in IS p.1517 Table3:
!          Pressure, specific internal energy , specific Cv, specific Cp, speed of sound,
!          specific entropy, specific helmholtz energy, specific gibbs energy
! @routine helmholtz_deriv.f90, helmholtz_dimless.f90
! @authors Yu Fang
! @date    10-02-2017
!--------------------------------------------------------------------------------------------------

MODULE properties
!
       USE def_constants 
!
       IMPLICIT NONE
! 
       PRIVATE
       PUBLIC :: pressure, inter_energy, heat_cap_v, heat_cap_p, &
&                sound_speed, entropy, helmho, gibbs,dp_dv !.....
       CONTAINS
!
!
!===============================================================================================
       SUBROUTINE pressure(T,v,p)  !Pa
!===============================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(pr) :: T,v
        REAL(pr) :: p              
!LOCAL
        REAL(pr) :: rho,delt
!
!for helmholtz_deriv
!
        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!Pre-compute
!
        rho      = 1_pr  / v
        delt     = 1_pr   / (rho_cr*v)
!
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
        p = (1_pr + delt * phi_r_ddelt) * R * T * rho
!
       END SUBROUTINE pressure
!
!=================================================================================================
       SUBROUTINE inter_energy(T,v,e) !J/kg
!=================================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(pr) :: T,v,e
!LOCAL
        REAL(pr) :: tau
!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau
!
!Pre-compute

        tau      = T_cr  / T
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr
!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau
!
!        print*, T,v
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
!
        e = tau * ( phi0_dtau + phi_r_dtau) * R * T
!
!        print*, 'phi_r_dtau(properties) = ',phi_r_dtau
!
       END SUBROUTINE inter_energy
!
!============================================================================================================
       SUBROUTINE heat_cap_v(T,v,cv)  !J/kgK
!============================================================================================================
        IMPLICIT NONE
!
!IN/OUT
        REAL(pr) :: T,v,cv
!LOCAL
        REAL(pr) :: tau, tau2
!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau
!
!Pre-compute

        tau      = T_cr  / T
        tau2     = tau   * tau
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr


        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
!
!
        cv = -tau2 * ( phi0_ddtau + phi_r_ddtau ) * R

!       print*, phi0_ddtau , phi_r_ddtau , tau2
!
!
       END SUBROUTINE heat_cap_v
!
!================================================================================================================        
       SUBROUTINE heat_cap_p(T,v,cp)  !J/kgK
!================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,cp
!LOCAL
        REAL(pr) :: rho,delt,tau,tau2,delt2

!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau


!Pre-compute

        rho      = 1_pr  / v
        delt     = rho   / rho_cr
        tau      = T_cr  / T
        tau2     = tau   * tau
        delt2    = delt  * delt
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr

        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
!
        cp = ( -tau2   * ( phi0_ddtau + phi_r_ddtau ) + &
             (1_pr+delt*phi_r_ddelt-delt*tau*phi_r_ddeltdtau)**(2_pr) / &
             (1_pr+2_pr*delt*phi_r_ddelt+delt2*phi_r_dddelt) ) * R
!
!
       END SUBROUTINE heat_cap_p
!
!=============================================================================================================
       SUBROUTINE sound_speed(T,v,c)   !m/s
!=============================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,c
!LOCAL
        REAL(pr) :: rho,delt,tau,tau2,delt2,c_2

!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau


!Pre-compute

        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        tau      = T_cr  / T
        tau2     = tau   * tau
        delt2    = delt  * delt

!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr

        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau,
!       phi_r_ddeltdtau
!
!
        c_2 = ( 1_pr + 2_pr*delt*phi_r_ddelt + delt2*phi_r_dddelt - &
              (1_pr + delt*phi_r_ddelt  - delt*tau*phi_r_ddeltdtau)**(2_pr) / &
              ( tau2 * (phi0_ddtau+phi_r_ddtau) ) ) * R * T
!        print*, tau2 * (phi0_ddtau+phi_r_ddtau)
!         print*, delt*phi_r_dddelt,delt*tau* phi_r_ddeltdtau
!         print*, (1_pr + delt*phi_r_ddelt  - delt*tau*phi_r_ddeltdtau)
!         print*, delt*phi_r_ddelt
        c = abs( sqrt(c_2) )
!
!
       END SUBROUTINE sound_speed
!
!
!====================================================================================================================
       SUBROUTINE entropy(T,v,s)  !J/kgK
!====================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,s
!LOCAL
        REAL(pr) :: tau
!for helmholtz_deriv

        REAL(pr) :: phi_r_ddelt, phi_r_dtau,&
                    phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau

!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau


!Pre-compute

        tau      = T_cr  / T
!
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr
      phi_r_ddelt= 0_pr
      phi_r_dtau = 0_pr
      phi_r_dddelt = 0_pr
      phi_r_ddtau = 0_pr
      phi_r_ddeltdtau = 0_pr
 
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!       print*, phi_0,phi_r,phi0_ddelt,phi0_dtau,phi0_dddelt,phi0_ddtau

        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
                            phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!       print*, phi_r_ddelt, phi_r_dtau,phi_r_dddelt, phi_r_ddtau,
!       phi_r_ddeltdtau
!
!
        s = (tau * ( phi0_dtau + phi_r_dtau) - phi_0 - phi_r) * R
!
!
       END SUBROUTINE entropy
!
!================================================================================================================
       SUBROUTINE helmho(T,v,h_helmho) !J/kg
!================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,h_helmho
!LOCAL
        REAL(pr) :: phi
!for helmholtz_dimless

        REAL(pr) :: phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                    phi0_dddelt,phi0_ddtau
!
!initialisation
      phi_0 = 0_pr
      phi_r = 0_pr
      phi0_ddelt = 0_pr
      phi0_dtau = 0_pr
      phi0_dddelt = 0_pr
      phi0_ddtau = 0_pr

!
        CALL helmholtz_dimless (T,v,phi_0,phi_r,phi0_ddelt,phi0_dtau,&
                              phi0_dddelt,phi0_ddtau)
!
!
        phi = phi_0 + phi_r
! print*, 'propoerties phi0, phir',phi_0, phi_r
! phi is a dimensionless magnitude, R is in [J/kg/K], T in [K]:
!
        h_helmho   = phi * R * T                      ! [J/kg]
!
!
       END SUBROUTINE helmho
!
!
!=================================================================================================================
       SUBROUTINE gibbs(T,v,g) !J/kg
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,g
!LOCAL
        REAL(pr) :: h_helmho, p
!
!
        CALL helmho(T,v,h_helmho)
        CALL pressure(T,v,p)
!
!
        g = h_helmho + p * v
!
!
       END SUBROUTINE gibbs
!
!
!==================================================================================================================
      SUBROUTINE dp_dv(T,v,deriv) !Pa.kg/m3
!=================================================================================================================
        IMPLICIT NONE
!IN/OUT
        REAL(pr) :: T,v,deriv
!LOCAL
        REAL(pr) :: rho, delt, phi_r_ddelt, phi_r_dtau,&
&                   phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau
!
        rho      = 1_pr  / v
        delt     = 1_pr  / (rho_cr*v)
        phi_r_ddelt = 0_pr
        phi_r_dtau  = 0_pr
        phi_r_dddelt= 0_pr
        phi_r_ddtau = 0_pr
        phi_r_ddeltdtau = 0_pr


!
        CALL helmholtz_deriv( T, v, phi_r_ddelt, phi_r_dtau,&
&                             phi_r_dddelt, phi_r_ddtau, phi_r_ddeltdtau )
!        
        deriv = (-1_pr - 2_pr*delt*phi_r_ddelt - phi_r_dddelt*delt*delt) *R *T *rho *rho


      END SUBROUTINE dp_dv


END MODULE properties
