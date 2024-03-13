module mod_treat_boundary_2D
   use mod_kinds
   use mod_field_2D, only: field
   use mod_co_func_2D
   implicit none
   private
   public :: treat_Dirchlet_boundary
contains
   subroutine treat_Dirchlet_boundary(Dirichlet_boundary_function, &
                                      A, b, boundarynodes, M_basis)
      real(wp), external :: Dirichlet_boundary_function
      real(wp), intent(inout) :: A(:, :)
      real(wp), intent(inout) :: b(:)
      integer, intent(in) :: boundarynodes(:, :)
      real(wp) :: M_basis(:, :)

      integer :: num_of_boundarynodes, k_index, index

      num_of_boundarynodes = size(boundarynodes, 2)
      !! Alogrithm III
      do k_index = 1, num_of_boundarynodes

         if (boundarynodes(1, k_index) == -1) then
            index = boundarynodes(2, k_index)
            A(index, :) = 0
            A(index, index) = 1
            b(index) = Dirichlet_boundary_function(M_basis(1, index), M_basis(2, index))

         end if

      end do
   end subroutine treat_Dirchlet_boundary
end module mod_treat_boundary_2D
