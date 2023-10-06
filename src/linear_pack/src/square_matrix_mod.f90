module square_matrix_mod
   implicit none
   private

   integer, parameter, public :: dp = selected_real_kind(15, 307)

   type, abstract, public :: square_matrix
   contains
      procedure(getm), deferred :: get_matrix
      procedure(gets), deferred :: get_size
      procedure, private :: mat_mult
      generic :: operator(*) => mat_mult
      procedure(solve), deferred :: inv_mat_mult
   end type square_matrix

   abstract interface
      pure function getm(this) result(m)
         import square_matrix
         import dp
         class(square_matrix), intent(in) :: this
         real(dp), dimension(this%get_size(), this%get_size()) :: m
      end function getm

      pure function gets(this) result(s)
         import square_matrix
         class(square_matrix), intent(in) :: this
         integer :: s
      end function gets

      function solve(this, rhs) result(solution)
         ! TODO: Consider returning a derived type with information on
         ! error, etc.
         import square_matrix
         import dp
         class(square_matrix), intent(inout) :: this
         real(dp), dimension(:), intent(in) :: rhs
         real(dp), dimension(size(rhs)) :: solution
      end function solve
   end interface

contains

   pure function mat_mult(this, rhs) result(product)
      class(square_matrix), intent(in) :: this
      real(dp), dimension(:), intent(in) :: rhs
      real(dp), dimension(size(rhs)) :: product
      real(dp), allocatable, dimension(:, :) :: mat
      if (this%get_size() /= size(rhs)) error stop "Matrix and array of different sizes"
      product = matmul(this%get_matrix(), rhs)
   end function mat_mult

end module square_matrix_mod

