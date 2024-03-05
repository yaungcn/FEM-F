module mod_basis_func_2D
   use mod_field_2D, only: field
   use mod_local_basis_func_2D, only: local_basis_func_2D
   use quadrature_module, only: wp => quadrature_wp
   implicit none

   private
   public :: trial_basis_func, test_basis_func

            !! basis_type: the type of the FE.
            !!             201-2D linear FE.
            !!             202-2D quadratic FE.
            !! basis_index: the index of basis function to specify
            !!              which basis function we want to use.
            !! der_x,der_y: derivative degree, the derivative degree of
            !!              the FE basis function.

   abstract interface
      pure real(wp) function trial_test_basis_func(field_info, vertices, &
                                                   basis_index, x, y, der_x, der_y)
         import :: field, wp
         class(field), intent(in) :: field_info
         real(wp), intent(in) :: vertices(2, 3)
         integer, intent(in) :: basis_index
         real(wp), intent(in) :: x, y
         integer, intent(in) :: der_x, der_y
      end function trial_test_basis_func
   end interface

contains

   pure real(wp) function trial_basis_func(field_info, vertices, &
                                           basis_index, x, y, der_x, der_y)
      class(field), intent(in) :: field_info
      real(wp), intent(in) :: vertices(2, 3)
      integer, intent(in) :: basis_index
      real(wp), intent(in) :: x, y
      integer, intent(in) :: der_x, der_y

      trial_basis_func = local_basis_func_2D(x, y, vertices, &
                                             field_info%trial_basis_type, &
                                             basis_index, der_x, der_y)

   end function trial_basis_func

   pure real(wp) function test_basis_func(field_info, vertices, &
                                          basis_index, x, y, der_x, der_y)
      class(field), intent(in) :: field_info
      real(wp), intent(in) :: vertices(2, 3)
      integer, intent(in) :: basis_index
      real(wp), intent(in) :: x, y
      integer, intent(in) :: der_x, der_y

      test_basis_func = local_basis_func_2D(x, y, vertices, &
                                            field_info%test_basis_type, &
                                            basis_index, der_x, der_y)
   end function test_basis_func
end module mod_basis_func_2D
