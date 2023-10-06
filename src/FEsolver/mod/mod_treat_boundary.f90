module mod_treat_boundary
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private
   public :: treat_Dirchlet_boundary
contains
   subroutine treat_Dirchlet_boundary(A, b)
      real(wp), allocatable, intent(in) :: A
      real(wp), allocatable, intent(inout) :: b

   end subroutine treat_Dirchlet_boundary
end module mod_treat_boundary
