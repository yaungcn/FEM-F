module mod_assemble_2D
   use mod_co_func_2D
   use mod_local_basis_func_2D
   use mod_gauss_quad_2D
   use mod_field_2D, only: field
   use quadrature_module, wp => quadrature_wp
   implicit none
   private
   public :: assemble_matirx_2D, assemble_vector_2D

   abstract interface
      pure function local_basis_func(x, y, vertices, basis_type, basis_index, der_x, der_y) result(res)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp), intent(in) :: vertices(2, 3)
         integer, intent(in) :: basis_type, basis_index, der_x, der_y
         real(wp) :: res
      end function local_basis_func
   !! add pure
      pure function cofunc_1d(x)
         import :: wp
         real(wp), intent(in) :: x
         real(wp) :: cofunc_1d
      end function cofunc_1d
   !! add pure, 2D
      pure function cofunc_2d(x, y)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp) :: cofunc_2d
      end function cofunc_2d
   !! add pure, 2D, time
      pure function cofunc_2d_t(x, y, t)
         import :: wp
         real(wp), intent(in) :: x, y, t
         real(wp) :: cofunc_2d_t
      end function cofunc_2d_t
   end interface

contains
   !! triangle_fe
   subroutine assemble_matirx_2D(A, co_func, lb_func, &
                                 M_partition, T_partition, &
                                 T_basis_trial, T_basis_test, &
                                 der_x_trial, der_y_trial, &
                                 der_x_test, der_y_test, &
                                 field_info)
      !> Alogrithm I-3
      real(wp), intent(inout) :: A(:, :)
      procedure(cofunc_2d), pointer, intent(in) :: co_func
                              !! co_func: The coefficient function
      procedure(local_basis_func), pointer, intent(in) :: lb_func
      real(wp), intent(in)    :: M_partition(:, :)
      integer, intent(in)     :: T_partition(:, :), &
                                 T_basis_trial(:, :), &
                                 T_basis_test(:, :)
      integer, intent(in)     :: der_x_trial, der_y_trial, &
                                 der_x_test, der_y_test
      type(field), intent(in) :: field_info

      real(wp) :: vertices(2, 3)
      !! In 2D, vertices is array with 2\times3 value.
      type(func_to_quad_2d) :: quadfunc
      real(wp) :: ans
      integer :: ierr
      real(wp) :: err

      integer  :: index_n, index_trial, index_test

      do index_n = 1, field_info%number_of_elements
         vertices(:, :) = M_partition(:, T_partition(:, index_n))
         do index_trial = 1, field_info%number_of_local_basis_trial
            do index_test = 1, field_info%number_of_local_basis_test
               call quadfunc%init_funcs(co_func, lb_func, &
                                        field_info, vertices, &
                                        der_x_trial, der_y_trial, index_trial, &
                                        der_x_test, der_y_test, index_test)
               call quadfunc%initialize(afunc, &
                                        vertices(1, 1), vertices(1, 2), &
                                        vertices(2, 1), vertices(2, 3), &
                                        field_info%tolerance, &
                                        field_info%tolerance, &
                                        field_info%Gauss_point_number, &
                                        field_info%Gauss_point_number)
               call quadfunc%integrate(ans, ierr, err)
               A(T_basis_test(index_test, index_n), &
                 T_basis_trial(index_trial, index_n)) &
                  = A(T_basis_test(index_test, index_n), T_basis_trial(index_trial, index_n)) + ans
            end do
         end do
      end do
   end subroutine assemble_matirx_2D

   subroutine assemble_vector_2D(b, co_func, lb_func, &
                                 M_partition, T_partition, &
                                 T_basis_test, &
                                 der_x_test, der_y_test, &
                                 field_info)
      !> Alogrithm II-3
      real(wp), intent(inout) :: b(:, :)
      procedure(cofunc_2d), pointer, intent(in) :: co_func
                              !! co_func: The coefficient function
      procedure(local_basis_func), pointer, intent(in) :: lb_func
      real(wp), intent(in)    :: M_partition(:, :)
      integer, intent(in)     :: T_partition(:, :), &
                                 T_basis_test(:, :)
      integer, intent(in)     :: der_x_test, der_y_test
      type(field), intent(in) :: field_info

      real(wp) :: vertices(2, 3)
      !! In 1D, vertices is array with 2 value.
      type(func_to_quad_2d) :: quadfunc
      real(wp) :: ans
      integer :: ierr
      real(wp) :: err

      integer  :: index_n, index_test

      do index_n = 1, field_info%number_of_elements
         vertices(:, :) = M_partition(:, T_partition(:, index_n))

         do index_test = 1, field_info%number_of_local_basis_test
            call quadfunc%init_funcs(co_func, lb_func, &
                                     field_info, vertices, &
                                     der_x_test, der_y_test, index_test, &
                                     der_x_test, der_y_test, index_test)
            call quadfunc%initialize(bfunc, &
                                     vertices(1, 1), vertices(1, 2), &
                                     vertices(2, 1), vertices(2, 3), &
                                     field_info%tolerance, &
                                     field_info%tolerance, &
                                     field_info%Gauss_point_number, &
                                     field_info%Gauss_point_number)
            call quadfunc%integrate(ans, ierr, err)
            b(T_basis_test(index_test, index_n), 1) = &
               b(T_basis_test(index_test, index_n), 1) + ans
         end do

      end do
   end subroutine assemble_vector_2D
end module mod_assemble_2D
