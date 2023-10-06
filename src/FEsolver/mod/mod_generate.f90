module mod_generate
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private
   public :: generate_info_matrix, generate_boundarynodes
contains
   subroutine generate_info_matrix(P, T)
      real(wp), allocatable, intent(inout) :: P, T
   end subroutine generate_info_matrix

   subroutine generate_boundarynodes(boundarynodes)
      integer, allocatable, intent(inout) :: boundarynodes

   end subroutine generate_boundarynodes
end module mod_generate
