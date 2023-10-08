module mod_assemble
   use mod_field, only: field
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private
   public :: assemble_matirx_1D, assemble_vector_1D
contains
   subroutine assemble_matirx_1D(A, field_info)
      real(wp), allocatable, intent(inout) :: A
      type(field), intent(in) :: field_info
   end subroutine assemble_matirx_1D

   subroutine assemble_vector_1D(b, field_info)
      real(wp), allocatable, intent(inout) :: b
      type(field), intent(in) :: field_info
   end subroutine assemble_vector_1D
end module mod_assemble
