module mod_linear_solver
   use quadrature_module, only: wp => quadrature_wp
   use linear_pack
   implicit none
   private
   public :: solver
contains
   subroutine solver(A, b, solution)
      real(wp), intent(in) :: A(:, :), b(:, :)
      real(wp), intent(inout) :: solution(:, :)
      type(general_square_matrix) :: linear_solver
      real(wp), dimension(:), allocatable :: b_solver
      real(wp), dimension(:), allocatable :: x_result
      integer :: size_n

      size_n = size(A, 1)

      allocate (b_solver(size_n), x_result(size_n))
      b_solver(:) = b(:, 1)
      linear_solver = general_square_matrix(A)
      x_result = linear_solver%inv_mat_mult(b_solver)

      solution(:, 1) = x_result(:)
   end subroutine solver
end module mod_linear_solver
