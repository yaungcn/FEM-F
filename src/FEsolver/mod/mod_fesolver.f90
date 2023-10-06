module mod_fesolver
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private
   public :: solver
contains
   subroutine solver(A, b, solution)
      real(wp), intent(in) :: A, b
      real(wp), allocatable, intent(inout) :: solution
   end subroutine solver
end module mod_fesolver
