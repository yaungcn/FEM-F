module mod_generate_1D
   use mod_field_1D, only: field
   use mod_kinds
   implicit none
   private
   public :: generate_info_matrix, generate_boundarynodes
contains
   subroutine generate_info_matrix(M, T, field_info)
      real(wp), allocatable, intent(inout) :: M(:, :)
      integer, allocatable, intent(inout) :: T(:, :)
      type(field), intent(in) :: field_info
      integer :: N, index

      if (field_info%basis_type == 101) then
         N = int((field_info%right - field_info%left)/field_info%h_partition)
         allocate (M(1, N + 1), T(2, N))

         do index = 1, N + 1
            M(1, index) = field_info%left + (index - 1)*field_info%h_partition
         end do

         do index = 1, N
            T(1, index) = index
            T(2, index) = index + 1
         end do
      end if
   end subroutine generate_info_matrix

   function generate_boundarynodes(N_basis)
      !! TODO: generater the information matrix for boundary nodes.
      !! Needs to be modified for different boundary conditions.
      !! Needs to be modified if we change te index of nodes and
      !!       elements in "generate_info_matrix".
      integer, allocatable :: generate_boundarynodes(:, :)
      integer, intent(in) :: N_basis

      allocate (generate_boundarynodes(3, 2))
      !! boundarynodes(1,k): specify the type of the kth boundary node.
      !!                   =-1: Dirichlet boundary node;
      !!                   =-2: Neumann boundary node;
      !!                   =-3 Robin boundary node.
      !! boundarynodes(2,k): global index of the kth boundary node among all nodes of FE.
      !!                      That is, the index of FE is used here.
      !! boundarynodes(3,k): The normal direction of the kth boundary nodes.
      generate_boundarynodes(1, 1) = -1
      generate_boundarynodes(2, 1) = 1
      generate_boundarynodes(3, 1) = -1
      generate_boundarynodes(1, 2) = -1
      generate_boundarynodes(2, 2) = N_basis + 1
      generate_boundarynodes(3, 2) = 1

   end function generate_boundarynodes
end module mod_generate_1D
