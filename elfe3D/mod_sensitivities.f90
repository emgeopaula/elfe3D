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
  use sparse_matrix_operations

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for obtaining Jacobian J times vector and 
  !> transposed Jacobian JT times vector in COO format
  !---------------------------------------------------------------------
  subroutine compute_Jvec_JTvec(E, num_rec, freq, &
                                forward_data, inv_model, free_M_indices, &
                                dAdrho, dAdrhorow, dAdrhocol, &
                                Jrows, Jcols, Jvec, JTvec)

  ! INPUT
  integer, intent(in) :: E, num_rec
  real(kind=dp), dimension(:), intent(in) :: freq 
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

   ! allocation


   ! initialise


   ! compute pseudo forward problems

   ! compute Jvec and JTvec

     ! loop over frequencies

       ! loop over free_M

       ! to include model parameter transformation log10 in free element loop
       ! drhodm = (10.0_dp**inv_model(i_free_M)) * log(10.0_dp)

       ! compute dAdrhoE or do before?

       ! loop over edges

         ! caclulate JTvec (sum over frequencies and edges)

         ! loop over data entries/2 to calculate Jvec



  

      
  ! -------------------------------------------------------------------
  end subroutine compute_Jvec_JTvec

  !-------------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating the RHS for and the solution of 
  !> the pseudo forward problems, that have to be solved to
  !> calculate the gradients of the objective function
  !> (see Newman & Alumbaugh, 2000)
  !> original version in emilia (in subroutine Compute_1v2v,
  !> implemented by Paula Rulff, 2023)
  !-------------------------------------------------------------------------
  subroutine compute_pseudo_fwd(E, freq, &
                                recx, recy, recz, &
                                rec_el, el2ed, &
                                a_start, a_end, b_start, b_end, &
                                c_start, c_end, d_start, d_end,&
                                el2edl, ed_sign, Ve, mu, &
                                pseudo_v, pseudo_u)

    ! INPUT
    integer, intent(in) :: E
    real(kind=dp), dimension(:), intent(in) :: freq
    real(kind=dp), dimension(:), intent(in) :: recx, recy, recz
    integer, dimension(:), intent(in) :: rec_el
    integer, dimension(:,:), intent(in) :: el2ed
    real(kind=dp), dimension(:,:), intent(in) :: a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end
    real(kind=dp), dimension(:,:), intent(in) :: el2edl
    real(kind=dp), dimension(:,:), intent(in) :: ed_sign
    real(kind=dp), dimension(:), intent(in) :: Ve
    real(kind=dp), dimension(:), intent(in) :: mu



    ! OUTPUT

    ! additional fwd solutions of pseudo forward problem for JTvec
    complex(kind=dp), dimension(:,:), intent(inout) :: pseudo_v
    ! additional fwdwd solutions of pseudo forward problem for Jvec
    complex(kind=dp), dimension(:,:,:), intent(inout) :: pseudo_u 



    ! local data
    ! counters
    integer :: ifreq, irec, idxf, l
    ! allocation
    integer :: allo_stat
    ! edges of the element containing receiver
    integer, dimension(6) :: rec1_ed 
    ! factor 1/(iwmu)
    complex(kind=dp) :: factor_mag   
    ! factor for g_datum calculations
    complex(kind=dp) :: factor 
    ! interpolator vectors for E and H field
    complex(kind=dp), allocatable, dimension(:,:) :: eN, hN 
    real(kind=dp), dimension(3) :: grad_Lstart
    real(kind=dp), dimension(3) :: grad_Lend

    ! RHS per datum for both source polarisations, solution v per frequency
    complex(kind=dp), allocatable, dimension(:,:) :: g_datum!, v_freq, t, u_freq
    ! RHS for both source polarisations and all sources
    !complex(kind=dp), allocatable, dimension(:,:,:) :: g 

    !-------------------------------------------------------------------------------
    ! PR: check + and - signs because convention is different from emilia!
    ! allocation
    allocate (eN(E,3), &
              hN(E,3), &
              g_datum(E,2), &
              stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
          "compute_pseudo_fwd: error allocating arrays for interpolator vectors!")

    ! initialise
    ifreq = 1
    grad_Lstart = 0.0_dp
    grad_Lend = 0.0_dp
    eN = cmplx(0.0_dp, 0.0_dp)
    hN = cmplx(0.0_dp, 0.0_dp)
    factor = cmplx(0.0_dp, 0.0_dp)
    factor_mag = cmplx(0.0_dp, 0.0_dp)
    g_datum = cmplx(0.0_dp, 0.0_dp)

    ! data item counter
    idxf = 0

    do ifreq = 1,size(freq)

      !if (allocated(g)) g = cmplx(0.0_dp, 0.0_dp)
      ! initialise
      idxf = 0
      irec = 0

      do irec = 1, size(rec_el)

        ! calculate interpolator function for the components of the 
        ! electric and magnetic field at the current station

        ! find edges of current receiver earth element
        rec1_ed = el2ed(rec_el(irec),:)

        ! calculate factor for magnetic interpolator vector
        ! (1/-(i*w*mu))
        !!!! PR: check sign!!!
        factor_mag = cmplx(0.0_dp, 0.0_dp)
        factor_mag = cmplx(0.0_dp, &
                           -(1.0_dp/(2.0_dp*pi*freq(ifreq) &
                            *(mu(rec_el(irec))))), &
                           kind=dp)

        ! calculate interpolation functions for receiver element
        ! re-initialise
        grad_Lstart = 0.0_dp
        grad_Lend = 0.0_dp
        eN = cmplx(0.0_dp, 0.0_dp)
        hN = cmplx(0.0_dp, 0.0_dp)

        do l = 1,6 ! edge loop
           ! calculate grad Lstart and grad Lend vectors
           grad_Lstart = (/ b_start(rec_el(irec),l), &
                            c_start(rec_el(irec),l), &
                            d_start(rec_el(irec),l) /)

           grad_Lend = (/ b_end(rec_el(irec),l), &
                          c_end(rec_el(irec),l), &
                          d_end(rec_el(irec),l) /)

           ! calculate Nedelec basis function for edge l
           eN(rec1_ed(l),:) = cmplx(((((1.0_dp/(6.0_dp*Ve(rec_el(irec))))**2.0_dp)* &
                                      (a_start(rec_el(irec),l) &
                                     + b_start(rec_el(irec),l) &
                                     * recx(irec) &
                                     + c_start(rec_el(irec),l) &
                                     * recy(irec) &
                                     + d_start(rec_el(irec),l) &
                                     * recz(irec))*grad_Lend &
                
                                     - &
          
                                      ((1.0_dp/(6.0_dp*Ve(rec_el(irec))))**2.0_dp)* &
                                       (a_end(rec_el(irec),l) &
                                      + b_end(rec_el(irec),l) &
                                      * recx(irec) &
                                      + c_end(rec_el(irec),l) &
                                      * recy(irec) &
                                      + d_end(rec_el(irec),l) &
                                      * recz(irec))*grad_Lstart) &
                
                                     * el2edl(rec_el(irec),l) &

                                     * ed_sign(rec_el(irec),l)), &
                                    kind=dp)


           hN(rec1_ed(l),:) = factor_mag * cmplx((((2.0_dp*el2edl(rec_el(irec),l)) &
                                                    /(6.0_dp*Ve(rec_el(irec)))**2.0_dp)* &
                                                 (/(c_start(rec_el(irec),l) &
                                                  * d_end(rec_el(irec),l) &
                                                  - d_start(rec_el(irec),l) &
                                                  * c_end(rec_el(irec),l)), &
                                                   (d_start(rec_el(irec),l) &
                                                  * b_end(rec_el(irec),l)&
                                                  - b_start(rec_el(irec),l) &
                                                  * d_end(rec_el(irec),l)), &
                                                   (b_start(rec_el(irec),l) &
                                                  * c_end(rec_el(irec),l) &
                                                  - c_start(rec_el(irec),l) &
                                                  * b_end(rec_el(irec),l))/) &
                                                  * ed_sign(rec_el(irec),l)), &
                                                  kind=dp)
        end do



      end do ! receiver loop

    end do ! frequency loop


    ! deallocate local arrays
    if (allocated(eN)) deallocate(eN)
    if (allocated(hN)) deallocate(hN)
    !if (allocated(g)) deallocate(g)
    if (allocated(g_datum)) deallocate(g_datum)
    !if (allocated(v_freq)) deallocate(v_freq)
    !if (allocated(t)) deallocate(t)
    !if (allocated(u_freq)) deallocate(u_freq)




    end subroutine compute_pseudo_fwd
   

end module sensitivities