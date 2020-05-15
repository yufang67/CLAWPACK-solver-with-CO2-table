!===========================================================================
!@ Module to create the grid in (e,v) domaine for LL, LH, R, HT and TP 
!@ The values of the properties (e,v,T,P,c) are computed by EoS Span-Wagner
!@ and available at each node of the grid.
!==========================================================================
MODULE Grid
!
!        USE def_constants
!        USE def_variables
!        USE non_linear_solvers
!      
        IMPLICIT NONE
!
!        PRIVATE
!
!        PUBLIC :: MAKE_GRID !
CONTAINS
!       
!
!
!==========================================================================
        SUBROUTINE MAKE_GRID()
!==========================================================================
!                
                CALL grid_construction_left_low
                CALL grid_construction_left_high
                CALL grid_construction_right
                CALL grid_construction_high_temperature
                CALL saturation_curve
                CALL grid_construction_TPL
                CALL grid_construction_TPM
                CALL grid_construction_TPH
!
        END SUBROUTINE MAKE_GRID
!
!
END MODULE Grid
