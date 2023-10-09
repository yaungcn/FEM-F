module mod_fesolver
   use quadrature_module, only: wp => quadrature_wp
   use linear_pack
   implicit none
   private
   public :: solver
contains
   subroutine solver(A, b, solution)
      real(wp), intent(in) :: A(:, :), b(:, :)
      real(wp), allocatable, intent(inout) :: solution(:, :)
      solution = A/b
   end subroutine solver
end module mod_fesolver
