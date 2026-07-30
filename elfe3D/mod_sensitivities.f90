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
  use solvers
  use sparse_matrix_operations

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for obtaining Jacobian J times vector and 
  !> transposed Jacobian JT times vector in COO format
  !---------------------------------------------------------------------
  subroutine compute_Jvec_JTvec(E, num_rec, freq, rec_el, &
                                forward_data, inv_model, free_M_indices, &
                                dAdrho, dAdrhorow, dAdrhocol, &
                                primal_solution, &
                                pseudo_v, pseudo_u, &
                                Jrows, Jcols, Jvec, JTvec)

  ! INPUT
  integer, intent(in) :: E, num_rec
  real(kind=dp), dimension(:), intent(in) :: freq 
  integer, dimension(:), intent(in) :: rec_el
  real(kind=dp), dimension(:), intent(in)  :: forward_data, inv_model
  integer, dimension(:), intent(in)  :: free_M_indices
  complex(kind=dp), dimension(:,:), intent(in) :: dAdrho
  integer, dimension(:,:), intent(in) :: dAdrhorow, dAdrhocol
  complex(kind=dp), dimension(:,:), intent(in) :: primal_solution
  ! additional fwd solutions of pseudo forward problem for JTvec
  complex(kind=dp), dimension(:,:), intent(in) :: pseudo_v
  ! additional fwd solutions of pseudo forward problem for Jvec
  complex(kind=dp), dimension(:,:,:), intent(in) :: pseudo_u 


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
  ! local product of dAdrho * solution
  complex(kind=dp), allocatable, dimension(:) :: dAdrhoE
  complex(kind=dp) :: dZdm
  real(kind=dp) :: drhodm 
  real(kind=dp) :: w
  integer :: i_free_M, iedge, ifreq, irec, i
  integer :: allo_stat
     
  !-------------------------------------------------------------------
   ! matrix entry indices
   Jrows = [(i, i = 1, size(forward_data))]
   Jcols = [(i, i = 1, size(inv_model))] ! PR (free_M_indices?, or 1-size inv_model?)

   ! allocation
   allocate (dAdrhoE(E), stat = allo_stat)
   call allocheck(log_unit, allo_stat, &
                  "Compute_JvecJTvec: error allocating dAdrhoE array!")


   ! initialise
   dAdrhoE = cmplx(0.0_dp, 0.0_dp, kind=dp)
   w = 0.0_dp
   drhodm = 0.0_dp


   ! compute Jvec and JTvec

   ! loop over frequencies
   do ifreq = 1,size(freq)

     ! define angular frequency
     w = 2.0_dp*pi*freq(ifreq)

     ! loop over free_M
     do i_free_M = 1, size(inv_model)

       ! to include model parameter transformation log10 in free element loop
       drhodm = (10.0_dp**inv_model(i_free_M)) * log(10.0_dp)

       ! compute dAdrho*Solution (including current angular frequency)
       dAdrhoE(:) = COMPSPARSEMUL(E,36,dAdrho(i_free_M,:)*w, &
                                       dAdrhorow(i_free_M,:), &
                                       dAdrhocol(i_free_M,:), &
                                       primal_solution(ifreq,:))

       ! loop over edges
       do iedge = 1, E

         ! caclulate JTvec (sum over frequencies and edges)
         JTvec(i_free_M) = JTvec(i_free_M) &
                           + real(pseudo_v(ifreq,iedge) * dAdrhoE(iedge))

         ! loop over data entries/2 to calculate Jvec
         ! PR: update for more than one data component
         do irec= 1, size(rec_el)

           ! caclulate "dZdm"
           dZdm = (pseudo_u(ifreq,iedge,irec) &
                   * dAdrhoE(iedge)) &
                   * cmplx(drhodm, 0.0_dp, kind=dp)

           ! sum over edges and elements, split real and imag
           ! Real part of Jvec
           Jvec(irec) = Jvec(irec) &
                        + real(dZdm) * inv_model(i_free_M)
           ! Imaginary part of Jvec
           Jvec(size(rec_el)+irec) = Jvec(size(rec_el)+irec) &
                                     + aimag(dZdm) * inv_model(i_free_M)


         end do ! end do receiver loop

       end do ! free element loop

     end do ! free element loop

   end do ! frequency loop



   ! deallocate local arrays
   deallocate(dAdrhoE)
  

      
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
                                system_matrix, jsystem_matrix, &
                                isystem_matrix, &
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
    complex(kind=dp), dimension(:,:), intent(in) :: system_matrix
    integer, dimension(:,:), intent(in) :: jsystem_matrix, isystem_matrix



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

    ! arrays for pseudo forward solutions
    complex(kind=dp), allocatable, dimension(:) :: g_datum, g, pseudo_v_freq
    complex(kind=dp), allocatable, dimension(:,:) :: t, pseudo_u_freq, qq_array

    !-------------------------------------------------------------------------------
    ! PR: check + and - signs because convention is different from emilia!
    ! allocation
    allocate (eN(E,3), &
              hN(E,3), &
              stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
          "compute_pseudo_fwd: error allocating arrays for interp vectors!")

    allocate (g_datum(E), &
              g(E), &
              pseudo_v_freq(E), &
              stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
     "compute_pseudo_fwd: error allocating arrays for pseudo fwd computation!")

    !PR: adapt with more data per receiver
    allocate (t(E, size(rec_el)), &
              qq_array(E, size(rec_el)), &
              pseudo_u_freq(E, size(rec_el)), &
              stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
     "compute_pseudo_fwd: error allocating arrays for pseudo fwd computation 2!")

    ! initialise
    ifreq = 1
    grad_Lstart = 0.0_dp
    grad_Lend = 0.0_dp
    eN = cmplx(0.0_dp, 0.0_dp, kind=dp)
    hN = cmplx(0.0_dp, 0.0_dp, kind=dp)
    factor = cmplx(0.0_dp, 0.0_dp, kind=dp)
    factor_mag = cmplx(0.0_dp, 0.0_dp, kind=dp)
    g_datum = cmplx(0.0_dp, 0.0_dp, kind=dp)

    g = cmplx(0.0_dp, 0.0_dp, kind=dp)
    pseudo_v_freq = cmplx(0.0_dp, 0.0_dp, kind=dp)
    t = cmplx(0.0_dp, 0.0_dp, kind=dp)
    qq_array = cmplx(0.0_dp, 0.0_dp, kind=dp)
    pseudo_u_freq = cmplx(0.0_dp, 0.0_dp, kind=dp)


    ! data item counter
    ! PR use to have more data items per receiver than just Ex
    idxf = 0

    do ifreq = 1,size(freq)
      
      ! initialise
      if (allocated(g)) g = cmplx(0.0_dp, 0.0_dp, kind=dp)
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
        factor_mag = cmplx(0.0_dp, 0.0_dp, kind=dp)
        factor_mag = cmplx(0.0_dp, &
                           -(1.0_dp/(2.0_dp*pi*freq(ifreq) &
                            *(mu(rec_el(irec))))), &
                           kind=dp)

        ! calculate interpolation functions for receiver element
        ! re-initialise
        grad_Lstart = 0.0_dp
        grad_Lend = 0.0_dp
        eN = cmplx(0.0_dp, 0.0_dp, kind=dp)
        hN = cmplx(0.0_dp, 0.0_dp, kind=dp)

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
        end do ! edge loop

        ! PR: adapt to more than Ex component later, check sign
        g_datum = -eN(:,1)

        ! building vector for JTvec calculation
        ! for sensitivity test with perturbation method:
        qq_array(ifreq,irec) = (1.0_dp, 0.0_dp)

        ! sum up for all data at one frequency for JTvec calculation only
        ! PR: now only for Ex component!
        g(:) = g(:) + CONJG(qq_array(ifreq,irec)) * g_datum

        ! do not sum up data items for Jvec calculation
        ! PR: change irec to idata later and add other components
        t(:,irec) = g_datum

      end do ! receiver loop

      ! solve pseudo forward problem for JTvec, solution: pseudo_v
      pseudo_v_freq = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call mumps_solving(system_matrix(ifreq,:), &
                         jsystem_matrix(ifreq,:), &
                         isystem_matrix(ifreq,:), &
                         g(:), pseudo_v_freq, 1)

      pseudo_v(ifreq,:) = pseudo_v_freq
      call Write_Message (log_unit, '*******************************************')
      call Write_Message (log_unit, 'pseudo forward problem for JTvec solved')
      call Write_Message (log_unit, '*******************************************')

      ! solve pseudo forward problem for Jvec, solution: pseudo_u
      pseudo_u_freq = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call mumps_solving_multiple_RHS (system_matrix(ifreq,:), &
                                       jsystem_matrix(ifreq,:), &
                                       isystem_matrix(ifreq,:), &
                                       t, pseudo_u_freq, size(rec_el))

      call Write_Message (log_unit, '*******************************************')
      call Write_Message (log_unit, 'pseudo forward problems for Jvec solved')
      call Write_Message (log_unit, '*******************************************')

      pseudo_u(ifreq,:,:) = pseudo_u_freq

    end do ! frequency loop


    ! deallocate local arrays
    if (allocated(eN)) deallocate(eN)
    if (allocated(hN)) deallocate(hN)
    if (allocated(g)) deallocate(g)
    if (allocated(g_datum)) deallocate(g_datum)
    if (allocated(pseudo_v_freq)) deallocate(pseudo_v_freq)
    if (allocated(t)) deallocate(t)
    if (allocated(pseudo_u_freq)) deallocate(pseudo_u_freq)

    contains

      ! subroutine to read in observed data and data errors from datafile
      ! located in /in and specified in elfe3D_input.txt
      subroutine read_observed_data
        ! open elfe3D_input.txt and check for line with 'input_data_file'

        ! read in filename (d_obs.txt)

        ! open d_obs.txt

        ! read data and data errors (same structure as qq_array (ifreq, irec))

        
      end subroutine read_observed_data
    end subroutine compute_pseudo_fwd
   

end module sensitivities