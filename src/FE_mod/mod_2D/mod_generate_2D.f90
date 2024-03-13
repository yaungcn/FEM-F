module mod_generate_2D
   use mod_kinds
   use mod_field_2D, only: field
   use mod_message, only: info_print
   implicit none
   private

   public :: generate_info_matrix, generate_boundarynodes, generate_boundaryedges

contains
   subroutine generate_info_matrix(field_info, M, T, basis_type, verbose)
      type(field), intent(in) :: field_info
      real(wp), allocatable, intent(inout) :: M(:, :)
      !! M(i, j): the coordinate of the j-th node in the field.
      integer, allocatable, intent(inout) :: T(:, :)
      !! T(i, j): the index of the nodes of the elements.
      integer, optional, intent(in) :: basis_type
      !! if trial and test basis function are different, generated different M and T.
      logical, optional, intent(in) :: verbose

      if (present(basis_type)) then
         call MT_matrix(M, T, field_info, basis_type)
      else
         call MT_matrix(M, T, field_info, field_info%trial_basis_type)
      end if

      if (present(verbose)) then
         if (verbose) then
            call info_print("The M and T information matrix generated.")
         end if
      end if
   end subroutine generate_info_matrix

   subroutine generate_boundarynodes(field_info, boundarynodes, verbose)
      type(field), intent(in) :: field_info
      integer, allocatable, intent(inout) :: boundarynodes(:, :)
      !! boundarynodes(1,k)=boundary_type: specify the type of the kth boundary node.
      !!                   =-1: Dirichlet boundary node;
      !!                   =-2: Neumann boundary node;
      !!                   =-3: Robin boundary node.
      !! boundarynodes(2,k): global index of the kth boundary node among all nodes of FE.
      !!                      That is, the index of FE is used here.
      !! boundarynodes(3,k): The normal direction of the kth boundary nodes.
      logical, optional, intent(in) :: verbose

      call boundary_nodes(boundarynodes, field_info)
      if (present(verbose)) then
         if (verbose) then
            call info_print("The boundarynodes information matrix generated.")
         end if
      end if

   end subroutine generate_boundarynodes

   subroutine generate_boundaryedges(field_info, boundaryedges, verbose)
      type(field), intent(in) :: field_info
      integer, allocatable, intent(inout) :: boundaryedges(:, :)
      !! boundaryedges(1,k)=boundary_type: specify the type of the kth boundary node.
      !!                   =-1: Dirichlet boundary node;
      !!                   =-2: Neumann boundary node;
      !!                   =-3: Robin boundary node.
      !! boundaryedges(2,k): index of the element which contains the kth boundary edge.
      !! boundaryedges(3:4,k): index of the two end points of the kth boundary edge among all grid points, not the nodes of FE.
      logical, optional, intent(in) :: verbose

      call boundary_edges(boundaryedges, field_info)

      if (present(verbose)) then
         if (verbose) then
            call info_print("The boundaryedges information matrix generated.")
         end if
      end if

   end subroutine generate_boundaryedges

   pure subroutine MT_matrix(M, T, field_info, basis_type)
      real(wp), allocatable, intent(inout) :: M(:, :)
      !! M(i, j): the coordinate of the j-th node in the field.
      integer, allocatable, intent(inout) :: T(:, :)
      !! T(i, j): the index of the nodes of the elements.
      type(field), intent(in) :: field_info
      integer, intent(in) :: basis_type
      integer :: Nh_p, Nv_p, Nh_b, Nv_b
      integer :: index, index_col, index_row, N_nodes, N_elements

      Nh_p = field_info%Nh_partition
      Nv_p = field_info%Nv_partition
      Nh_b = field_info%Nh_basis
      Nv_b = field_info%Nv_basis

      select case (basis_type)
       case (201)
         !> For 2D linear FE, M_partition,T_partition are the same as M_basis,T_basis.
         !! N_nodes: the number of nodes in the field.
         !! N_elements: the number of elements in the field,
         !! M: the coordinate of the nodes in the field.
         !! T: the index of the nodes of the elements.
         select case (field_info%mesh_type)

          case (503)!! 503: triangular mesh information matrix.
            N_nodes = (Nh_p + 1)*(Nv_p + 1)
            N_elements = Nh_p*Nv_p*2 !! 2 triangle elements in each rectangle.
            allocate (M(2, N_nodes), T(3, N_elements))

            !> generate the information matrix M and T.
            do index_col = 1, Nh_p + 1
               do index_row = 1, Nv_p + 1
                  !! index of the node in the field.
                  index = (index_col - 1)*(Nv_p + 1) + index_row
                  M(1, index) = field_info%left + (index_col - 1)*field_info%h_partition
                  M(2, index) = field_info%bottom + (index_row - 1)*field_info%v_partition
               end do
            end do

            do index_col = 1, Nh_p
               do index_row = 1, Nv_p
                  !! index of the element in the field.
                  index = 2*index_row + 2*(index_col - 1)*Nv_p
                  !! down triangle.
                  T(1, index - 1) = (Nv_p + 1)*(index_col - 1) + index_row
                  T(2, index - 1) = T(1, index - 1) + Nv_p + 1
                  T(3, index - 1) = T(1, index - 1) + 1
                  !! up triangle.
                  T(1, index) = T(3, index - 1)
                  T(2, index) = T(3, index - 1) + Nv_p
                  T(3, index) = T(2, index) + 1
               end do
            end do

          case (504)!! 504: quadrilateral mesh information matrix.
            N_nodes = (Nh_p + 1)*(Nv_p + 1)
            N_elements = Nh_p*Nv_p
            allocate (M(2, N_nodes), T(4, N_elements))

            !> generate the information matrix M and T.
            do index_col = 1, Nh_p + 1
               do index_row = 1, Nv_p + 1
                  !! index of the node in the field.
                  index = (index_col - 1)*(Nv_p + 1) + index_row
                  M(1, index) = field_info%left + (index_col - 1)*field_info%h_partition
                  M(2, index) = field_info%bottom + (index_row - 1)*field_info%v_partition
               end do
            end do

            do index_col = 1, Nh_p
               do index_row = 1, Nv_p
                  !! index of the element in the field.
                  index = index_row + (index_col - 1)*Nv_p
                  !! index of the node in the element.
                  T(1, index) = (Nv_p + 1)*(index_col - 1) + index_row
                  T(2, index) = T(1, index) + Nv_p + 1
                  T(3, index) = T(2, index) + 1
                  T(4, index) = T(1, index) + 1

               end do
            end do
          case default
            error stop "Error: The mesh type is not supported."
         end select

       case (202)
         !> For 2D Langrange quadratic FE, M_partition,T_partition are different from M_basis,T_basis.
         !! N_nodes: the number of nodes in the field.
         !! N_elements: the number of elements in the field,
         !! M: the coordinate of the nodes in the field.
         !! T: the index of the nodes of the elements.

         select case (field_info%mesh_type)
          case (503)!! 503: triangular mesh information matrix.
            N_elements = Nh_p*Nv_p*2 !! 2 triangle elements in each rectangle.
            N_nodes = (Nh_b + 1)*(Nv_b + 1)
            allocate (M(2, N_nodes), T(6, N_elements))

            !> generate the information matrix M and T.
            do index_col = 1, (Nh_b + 1)
               do index_row = 1, (Nv_b + 1)
                  !! index of the node in the field.
                  index = (index_col - 1)*(Nv_b + 1) + index_row
                  M(1, index) = field_info%left + (index_col - 1)*field_info%h_partition/2
                  M(2, index) = field_info%bottom + (index_row - 1)*field_info%v_partition/2
               end do
            end do

            do index_col = 1, Nh_p
               do index_row = 1, Nv_p
                  !! index of the element in the field.
                  index = 2*index_row + 2*(index_col - 1)*Nv_p
                  !! down triangle.
                  T(1, index - 1) = 2*(Nv_b + 1)*(index_col - 1) + 2*index_row - 1
                  T(2, index - 1) = T(1, index - 1) + 2*(Nv_b + 1)
                  T(3, index - 1) = T(1, index - 1) + 2
                  T(4, index - 1) = T(1, index - 1) + Nv_b + 1
                  T(5, index - 1) = T(4, index - 1) + 1
                  T(6, index - 1) = T(1, index - 1) + 1
                  !! up triangle.
                  T(1, index) = T(3, index - 1)
                  T(2, index) = T(2, index - 1)
                  T(3, index) = T(2, index) + 2
                  T(4, index) = T(4, index - 1) + 1
                  T(5, index) = T(2, index) + 1
                  T(6, index) = T(4, index) + 1
               end do
            end do

          case (504)!! 504: quadrilateral mesh information matrix.
            N_elements = Nh_p*Nv_p
            N_nodes = (Nh_b + 1)*(Nv_b + 1)
            allocate (M(2, N_nodes), T(9, N_elements))

            !> generate the information matrix M and T.
            do index_col = 1, (Nh_b + 1)
               do index_row = 1, (Nv_b + 1)
                  !! index of the node in the field.
                  index = (index_col - 1)*(Nv_b + 1) + index_row
                  M(1, index) = field_info%left + (index_col - 1)*field_info%h_partition/2
                  M(2, index) = field_info%bottom + (index_row - 1)*field_info%v_partition/2
               end do
            end do

            do index_col = 1, Nh_p
               do index_row = 1, Nv_p
                  !! index of the element in the field.
                  index = index_row + (index_col - 1)*Nv_p
                  !! rectangle.
                  T(1, index) = 2*(Nv_b + 1)*(index_col - 1) + 2*index_row - 1
                  T(2, index) = T(1, index) + 2*(Nv_b + 1)
                  T(3, index) = T(2, index) + 2
                  T(4, index) = T(1, index) + 2
                  T(5, index) = T(1, index) + Nv_b + 1
                  T(6, index) = T(2, index) + 1
                  T(7, index) = T(5, index) + 2
                  T(8, index) = T(1, index) + 1
                  T(9, index) = T(5, index) + 1
               end do
            end do

          case default
            error stop "Error: The mesh type is not supported."
         end select

      end select
   end subroutine MT_matrix

   pure subroutine boundary_nodes(boundarynodes, field_info)
      !! TODO: generater the information matrix for boundary nodes.
      !! Needs to be modified for different boundary conditions.
      !! Needs to be modified if we change te index of nodes and
      !!       elements in "generate_info_matrix".
      integer, allocatable, intent(inout) :: boundarynodes(:, :)
      type(field), intent(in) :: field_info
      integer :: Nh_b, Nv_b
      integer :: N_boundarynodes, boundary_type
      integer :: index

      Nh_b = field_info%Nh_basis
      Nv_b = field_info%Nv_basis
      N_boundarynodes = 2*(Nh_b + Nv_b)
      boundary_type = -1

      allocate (boundarynodes(2, N_boundarynodes))
      !! boundarynodes(1,k)=boundary_type: specify the type of the kth boundary node.
      !!                   =-1: Dirichlet boundary node;
      !!                   =-2: Neumann boundary node;
      !!                   =-3: Robin boundary node.
      !! boundarynodes(2,k): global index of the kth boundary node among all nodes of FE.
      !!                      That is, the index of FE is used here.
      !! boundarynodes(3,k): The normal direction of the kth boundary nodes.
      boundarynodes(1, :) = boundary_type

      !! bottom boundary.
      do index = 1, Nh_b
         boundarynodes(2, index) = (index - 1)*(Nv_b + 1) + 1
      end do

      !! right boundary.
      do index = Nh_b + 1, Nh_b + Nv_b
         boundarynodes(2, index) = Nh_b*(Nv_b + 1) + index - Nh_b
      end do

      !! top boundary.
      do index = Nh_b + Nv_b + 1, 2*Nh_b + Nv_b
         boundarynodes(2, index) = (2*Nh_b + Nv_b + 2 - index)*(Nv_b + 1)
      end do

      !! left boundary.
      do index = 2*Nh_b + Nv_b + 1, N_boundarynodes
         boundarynodes(2, index) = 2*Nh_b + 2*Nv_b + 2 - index
      end do

   end subroutine boundary_nodes

   pure subroutine boundary_edges(boundaryedges, field_info)
      !! TODO: generater the information matrix for boundary nodes.
      !! Needs to be modified for different boundary conditions.
      !! Needs to be modified if we change te index of nodes and
      !!       elements in "generate_info_matrix".
      integer, allocatable, intent(inout) :: boundaryedges(:, :)
      type(field), intent(in) :: field_info
      integer :: Nh_p, Nv_p
      integer :: N_boundaryedges, boundary_type
      integer :: index

      Nh_p = field_info%Nh_partition
      Nv_p = field_info%Nv_partition
      N_boundaryedges = 2*(Nh_p + Nv_p)
      boundary_type = -1

      allocate (boundaryedges(4, N_boundaryedges))
      !! boundaryedges(1,k)=boundary_type: specify the type of the kth boundary node.
      !!                   =-1: Dirichlet boundary node;
      !!                   =-2: Neumann boundary node;
      !!                   =-3: Robin boundary node.
      !! boundaryedges(2,k): index of the element which contains the kth boundary edge.
      !! boundaryedges(3:4,k): index of the two end points of the kth boundary edge among all grid points, not the nodes of FE.
      boundaryedges(1, :) = boundary_type

      !! bottom boundary.
      do index = 1, Nh_p
         boundaryedges(2, index) = (index - 1)*2*Nv_p + 1
         boundaryedges(3, index) = (index - 1)*(Nv_p + 1) + 1
         boundaryedges(4, index) = index*(Nv_p + 1) + 1
      end do

      !! right boundary.
      do index = Nh_p + 1, Nh_p + Nv_p
         boundaryedges(2, index) = (Nh_p - 1)*2*Nv_p + 2*(index - Nh_p)
         boundaryedges(3, index) = Nh_p*(Nv_p + 1) + index - Nh_p
         boundaryedges(4, index) = Nh_p*(Nv_p + 1) + index - Nh_p + 1
      end do

      !! top boundary.
      do index = Nh_p + Nv_p + 1, 2*Nh_p + Nv_p
         boundaryedges(2, index) = (2*Nh_p + Nv_p + 1 - index)*2*Nv_p
         boundaryedges(3, index) = (2*Nh_p + Nv_p + 2 - index)*(Nv_p + 1)
         boundaryedges(4, index) = (2*Nh_p + Nv_p + 1 - index)*(Nv_p + 1)
      end do

      !! left boundary.
      do index = 2*Nh_p + Nv_p + 1, N_boundaryedges
         boundaryedges(2, index) = 2*(2*Nh_p + 2*Nv_p + 1 - index) - 1
         boundaryedges(3, index) = 2*Nh_p + 2*Nv_p + 2 - index
         boundaryedges(4, index) = 2*Nh_p + 2*Nv_p + 1 - index
      end do

   end subroutine boundary_edges
end module mod_generate_2D
