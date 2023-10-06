module mod_assemble
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private
   public :: assemble_matirx_1D, assemble_vector_1D
contains
   subroutine assemble_matirx_1D(A)
      real(wp), allocatable, intent(inout) :: A
   end subroutine assemble_matirx_1D

   subroutine assemble_vector_1D(b)
      real(wp), allocatable, intent(inout) :: b
   end subroutine assemble_vector_1D
end module mod_assemble
