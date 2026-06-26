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
  subroutine compute_Jvec_JTvec(forward_data, inv_model, free_M_indices, &
                                dAdrho, dAdrhorow, dAdrhocol, &
                                Jrows, Jcols, Jvec, JTvec)

  ! INPUT
  real(kind=dp), dimension(:), intent(in)  :: forward_data, inv_model
  integer, dimension(:), intent(in)  :: free_M_indices
  complex(kind=dp), dimension(:,:), intent(in) :: dAdrho
  integer, dimension(:,:), intent(in) :: dAdrhorow, dAdrhocol

  ! OUTPUT
  ! sensitivity output: Jacobian times vectors in COO format
  ! number of Jrows corresponds to size(forward_data)
  ! number of columns corresponds to the free model size (num_free_M)
  integer(kind=dp), dimension(:), intent(inout) :: Jrows
  integer(kind=dp), dimension(:), intent(inout) :: Jcols
  ! J*inv_model
  real(kind=dp), dimension(:), intent(inout) :: Jvec 
  ! JT*forward_data = (forward_data*JT)T
  ! doublecheck if forward_data should be d = (dobs-dpred)/epsilon
  real(kind=dp), dimension(:), intent(inout) :: JTvec 

  ! LOCAL VARIABLES
  real(kind=dp) :: drhodm 
  integer :: i_free_M
     

  !-------------------------------------------------------------------
   ! dummy entries
   Jvec = 666.0_dp
   JTvec = 444.0_dp
   Jrows = 6
   Jcols = 4


   ! to include model parameter transformation log10 in free element loop
   ! drhodm = (10.0_dp**inv_model(i_free_M)) * log(10.0_dp)

  ! -------------------------------------------------------------------
  end subroutine compute_Jvec_JTvec
   

end module sensitivities