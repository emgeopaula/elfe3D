!> @brief
!> Module of elfe3d containing subroutines to calcullate sensitivities
!> NEW for elfe3DINV
!!
!> written by Paula Rulff, 5/5/2026
!!
!> Last change: May 2026
!!
!> Copyright (C) Paula Rulff 2026
!!
module sensitivities

  use mod_util

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for obtaining Jacobian J times vector and 
  !> transposed Jacobian JT times vector in COO format
  !---------------------------------------------------------------------
  subroutine compute_Jvec_JTvec(Jrows, Jcols, Jvec, JTvec)
  ! add as input: data?, forward solution
  ! rho, num_free_M, free_M_indices, &
  ! Adrho, dAdrhorow, dAdrhocol, &

  ! INPUT
  ! PR: update

  ! OUTPUT
  ! sensitivity output: Jacobian times vectors in COO format
  integer(kind=dp), dimension(:), intent(inout) :: Jrows
  integer(kind=dp), dimension(:), intent(inout) :: Jcols
  ! number of Jrows corresponds to 2x data-size (Re; Im)
  ! number of columns corresponds to the model size (num_free_M)
  real(kind=dp), dimension(:), intent(inout) :: Jvec
  real(kind=dp), dimension(:), intent(inout) :: JTvec
     

  !-------------------------------------------------------------------
   ! dummy entries
   Jvec = 666.0_dp
   JTvec = 444.0_dp
   Jrows = 6
   Jcols = 4

  ! -------------------------------------------------------------------
  end subroutine compute_Jvec_JTvec
   

end module sensitivities