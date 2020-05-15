MODULE var_const

  IMPLICIT NONE
  
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
    real(8),parameter :: e_const = 5.0d6 
    integer,parameter :: Ngost   = 2
    integer,parameter :: Nxmax   = 1000
    integer,parameter :: Nymax   = 1000
    integer,parameter :: N_1d    = 100000




!-----------------------------------------------------------------
!
!-----------------------------------------------------------------

    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: Tguess
    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: guessP
    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: guessP_rpn
    real(8),dimension(1-Ngost:Nxmax+Ngost,1-Ngost:Nymax+Ngost) :: c_BC
    real(8),dimension(1-Ngost:N_1d+Ngost) :: guessp_1d
    real(8),dimension(1-Ngost:Nymax+Ngost) :: guessp_BC
    real(8),dimension(1-Ngost:Nymax+Ngost) :: vguess_in
    real(8),dimension(1-Ngost:Nymax+Ngost) :: eguess_out
    integer :: i_bc
    integer :: flag_diagno
    integer :: flag_out
    integer :: i_flag
    real(8) :: e_guess
    real(8) :: stepT
    real(8) :: stepx

END MODULE var_const
