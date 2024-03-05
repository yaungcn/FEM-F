module mod_assemble_1D
   use mod_co_func_1D
   use mod_basis_func_1D
   use mod_gauss_quad_1D
   use mod_field_1D, only: field
   use quadrature_module, wp => quadrature_wp
   implicit none
   private
   public :: assemble_matirx_1D, assemble_vector_1D

contains
   subroutine assemble_matirx_1D(A, co_func, &
                                 M_partition, T_partition, &
                                 T_basis_trial, T_basis_test, &
                                 num_of_elements, &
                                 num_of_trial_local_basis, num_of_test_local_basis, &
                                 trial_basis_type, trial_derivate_degree, &
                                 test_basis_type, test_derivate_degree, &
                                 field_info)
      real(wp), intent(inout) :: A(:, :)
      class(func), intent(in) :: co_func
                              !! co_func: The coefficient function
      real(wp), intent(in)    :: M_partition(:, :)
      integer, intent(in)     :: T_partition(:, :), &
                                 T_basis_trial(:, :), &
                                 T_basis_test(:, :)
      integer, intent(in)     :: num_of_elements, &
                                 num_of_trial_local_basis, &
                                 num_of_test_local_basis, &
                                 trial_basis_type, trial_derivate_degree, &
                                 test_basis_type, test_derivate_degree
      type(field), intent(in) :: field_info

      type(local_basis_1D) :: trial_basis, test_basis
      real(wp) :: vertices(2)
      !! In 1D, vertices is array with 2 value.
      type(quad_1D) :: gauss_quad
      real(wp) :: ans
      integer :: ierr
      real(wp) :: err

      integer  :: index_n, index_trial, index_test

      do index_n = 1, num_of_elements
         vertices(1) = M_partition(1, T_partition(1, index_n))
         vertices(2) = M_partition(1, T_partition(2, index_n))

         do index_trial = 1, num_of_trial_local_basis
            do index_test = 1, num_of_test_local_basis
               call trial_basis%init(vertices, trial_basis_type, index_trial, trial_derivate_degree)
               call test_basis%init(vertices, test_basis_type, index_test, test_derivate_degree)
               call gauss_quad%init(co_func, trial_basis, test_basis)
               call gauss_quad%initialize(quad_A_1D, &
                                          vertices(1), vertices(2), &
                                          field_info%tolerance, &
                                          field_info%Gauss_point_number)
               call gauss_quad%integrate(ans, ierr, err)

               A(T_basis_test(index_test, index_n), T_basis_trial(index_trial, index_n)) = &
                  A(T_basis_test(index_test, index_n), T_basis_trial(index_trial, index_n)) + ans
            end do
         end do
      end do
   end subroutine assemble_matirx_1D

   subroutine assemble_vector_1D(b, co_func, &
                                      M_partition, T_partition, &
                                      T_basis_test, num_of_elements, &
                                      num_of_test_local_basis, &
                                      test_basis_type, test_derivate_degree, &
                                      field_info)
      real(wp), intent(inout) :: b(:, :)
      class(func), intent(in) :: co_func
                              !! co_func: The coefficient function
      real(wp), intent(in)    :: M_partition(:, :)
      integer, intent(in)     :: T_partition(:, :), &
                                 T_basis_test(:, :)
      integer, intent(in)     :: num_of_elements, &
                                 num_of_test_local_basis, &
                                 test_basis_type, test_derivate_degree
      type(field), intent(in) :: field_info

      type(local_basis_1D) :: test_basis
      real(wp) :: vertices(2)
      !! In 1D, vertices is array with 2 value.
      type(quad_1D) :: gauss_quad
      real(wp) :: ans
      integer :: ierr
      real(wp) :: err

      integer  :: index_n, index_test

      do index_n = 1, num_of_elements
         vertices(1) = M_partition(1, T_partition(1, index_n))
         vertices(2) = M_partition(1, T_partition(2, index_n))

         do index_test = 1, num_of_test_local_basis

            call test_basis%init(vertices, test_basis_type, index_test, test_derivate_degree)
            call gauss_quad%init(co_func, test_basis, test_basis)
            call gauss_quad%initialize(quad_b_1D, &
                                       vertices(1), vertices(2), &
                                       field_info%tolerance, &
                                       field_info%Gauss_point_number)
            call gauss_quad%integrate(ans, ierr, err)

            b(T_basis_test(index_test, index_n), 1) = &
               b(T_basis_test(index_test, index_n), 1) + ans
         end do

      end do
   end subroutine assemble_vector_1D
end module mod_assemble_1D
