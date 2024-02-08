module mod_generate_2D
   use mod_field_2D, only: field
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private
   public :: generate_info_matrix, generate_boundarynodes
contains
   subroutine generate_info_matrix(M, T, field_info)
      real(wp), allocatable, intent(inout) :: M(:, :)
      !! M(i, j): the coordinate of the j-th node in the field.
      integer, allocatable, intent(inout) :: T(:, :)
      !! T(i, j): the index of the nodes of the elements.
      type(field), intent(in) :: field_info
      integer :: index, index_col, index_row, N_nodes, N_elements

      if (field_info%mesh_type == 503) then
         !! 503: triangular mesh information matrix.
         N_nodes = (field_info%Nh + 1)*(field_info%Nv + 1)
         N_elements = 2*field_info%Nh*field_info%Nv
         allocate (M(2, N_nodes), T(3, N_elements))

         do index_col = 1, field_info%Nh + 1
            do index_row = 1, field_info%Nv + 1
               !! index of the node in the field.
               index = (index_col - 1)*(field_info%Nv + 1) + index_row
               M(1, index) = field_info%left + (index_col - 1)*field_info%h_partition
               M(2, index) = field_info%bottom + (index_row - 1)*field_info%v_partition
            end do
         end do

         do index_col = 1, field_info%Nh
            do index_row = 1, field_info%Nv
               !! index of the element in the field.
               index = 2*index_row + 2*(index_col - 1)*field_info%Nv
               !! down triangle.
               T(1, index - 1) = (field_info%Nv + 1)*(index_col - 1) + index_row
               T(2, index - 1) = T(1, index - 1) + field_info%Nv + 1
               T(3, index - 1) = T(1, index - 1) + 1
               !! up triangle.
               T(1, index) = T(3, index - 1)
               T(2, index) = T(3, index - 1) + field_info%Nv
               T(3, index) = T(2, index) + 1
            end do
         end do

      else if (field_info%mesh_type == 504) then
         !! 504: quadrilateral mesh information matrix.
         N_nodes = (field_info%Nh + 1)*(field_info%Nv + 1)
         N_elements = field_info%Nh*field_info%Nv
         allocate (M(2, N_nodes), T(4, N_elements))

         do index_col = 1, field_info%Nh + 1
            do index_row = 1, field_info%Nv + 1
               !! index of the node in the field.
               index = (index_col - 1)*(field_info%Nv + 1) + index_row
               M(1, index) = field_info%left + (index_col - 1)*field_info%h_partition
               M(2, index) = field_info%bottom + (index_row - 1)*field_info%v_partition
            end do
         end do

         do index_col = 1, field_info%Nh
            do index_row = 1, field_info%Nv
               !! index of the element in the field.
               index = index_row + (index_col - 1)*field_info%Nv
               !! index of the node in the element.
               T(1, index) = (field_info%Nv + 1)*(index_col - 1) + index_row
               T(2, index) = T(1, index) + field_info%Nv + 1
               T(3, index) = T(2, index) + 1
               T(4, index) = T(1, index) + 1

            end do
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
end module mod_generate_2D
