program main
   use mod_kinds
   use M_attr
   use mod_FE_2D

   implicit none
   interface
      pure function local_basis_func(x, y, local_vertices, lbasis_type, basis_index, der_x, der_y) result(res)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp), intent(in) :: local_vertices(2, 3)
         integer, intent(in) :: lbasis_type, basis_index, der_x, der_y
         real(wp) :: res
      end function local_basis_func
      !! 1D
      pure function cofunc_1d(x)
         import :: wp
         real(wp), intent(in) :: x
         real(wp) :: cofunc_1d
      end function cofunc_1d
      !! 2D
      pure function cofunc_2d(x, y)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp) :: cofunc_2d
      end function cofunc_2d
      !! 2D, time
      pure function cofunc_2d_t(x, y, time)
         import :: wp
         real(wp), intent(in) :: x, y, time
         real(wp) :: cofunc_2d_t
      end function cofunc_2d_t
   end interface

   !> input information of field
   type(field) :: field_info

!> Field parameter input.
   !>>TODO: read from toml file.
   real(wp) :: left = 0.0_wp, right = 10.0_wp, bottom = 0.0_wp, top = 10.0_wp
   !! The problem domain is [left,right]*[bottom,top].
   integer :: Nh_parition = 100, Nv_parition = 100
   !! The number of partition in horizontal and vertical direction.
   ! real(wp) :: h_partition, v_partition
   !! The step size of the partition.
   integer :: Nh_basis, Nv_basis
   !! The number of FE basis functions in horizontal and vertical direction.
   integer :: Gauss_point_number = 6
   !! The number of Gauss Quadrature points.
   real(wp) :: tolerance = 1.0e-3_wp
   !! The tolerance of the error.
   !> basis_type: the type of the FE.
   integer :: basis_type = 201
   !! 201: 2D linear
   !! 202: 2D quadratic
   integer :: mesh_type = 503
   !! 503: 2D triangular mesh;
   !! 504: 2D quadrilateral mesh; !!!NOT IMPLEMENTED YET!!!

   real(wp), allocatable :: M(:, :), M_basis(:, :)
   integer, allocatable :: T(:, :), T_basis(:, :)
   integer, allocatable :: boundarynodes(:, :)
   integer, allocatable :: boundaryedges(:, :)

   real(wp), allocatable :: A(:, :), A1(:, :), A2(:, :)
   real(wp), allocatable :: b(:), solution(:)
   ! real(wp) :: vet(2, 3)

   !> sparse matrix and solver
   type(sparse_COO):: sparse_A, sparseA1, sparseA2

!> co-Function input.
   procedure(local_basis_func), pointer :: lbfunc => null()
   procedure(cofunc_2d), pointer :: cofunc => null()
   procedure(cofunc_2d), pointer :: cofuncf => null()
   lbfunc => local_basis_func_2D
   cofunc => function_c
   cofuncf => function_f

!> Field information initialization.
   select case (basis_type)
    case (201)
      Nh_basis = Nh_parition
      Nv_basis = Nv_parition
    case (202)
      Nh_basis = 2*Nh_parition
      Nv_basis = 2*Nv_parition
   end select

   call field_info%init(left, right, bottom, top, &
      Nh_parition, Nv_parition, &
      Nh_basis, Nv_basis, &
      Gauss_point_number=Gauss_point_number, &
      tolerance=tolerance, &
      trial_basis_type=basis_type, &
      test_basis_type=basis_type, &
      mesh_type=mesh_type)

   ! call field_info%print('(A,2F6.3)')

   !> Allocate memory to A, b and solution. Set to zero.
   allocate (A(field_info%number_of_nodes_fe, field_info%number_of_nodes_fe))
   allocate (A1(field_info%number_of_nodes_fe, field_info%number_of_nodes_fe))
   allocate (A2(field_info%number_of_nodes_fe, field_info%number_of_nodes_fe))
   allocate (b(field_info%number_of_nodes_fe))
   allocate (solution(field_info%number_of_nodes_fe))
   A1 = 0
   A2 = 0
   A = 0
   b = 0
   solution = 0

   !> Generate the mesh and basis information matrix.
   call generate_info_matrix(field_info, M, T, &
      basis_type=basis_type, verbose=.true.)
   call generate_info_matrix(field_info, M_basis, T_basis, &
      basis_type=basis_type, verbose=.true.)
   !> Generate the boundary nodes and edges matrix.
   call generate_boundarynodes(field_info, boundarynodes, verbose=.true.)
   call generate_boundaryedges(field_info, boundaryedges, verbose=.true.)

   !! Assemble the matrix and vector.
   !! Triangular mesh.
   call assemble_matirx_2D(A1, cofunc, lbfunc, &
      M, T, &
      T_basis, T_basis, &
      1, 0, &
      1, 0, &
      field_info)
   call assemble_matirx_2D(A2, cofunc, lbfunc, &
      M, T, &
      T_basis, T_basis, &
      0, 1, &
      0, 1, &
      field_info)

   A = A1 + A2
   deallocate (A1, A2)
   call assemble_vector_2D(b, cofuncf, lbfunc, &
      M, T, &
      T_basis, &
      0, 0, &
      field_info)
   call treat_Dirchlet_boundary(function_g, A, b, boundarynodes, M)

   !! Solve the linear system.
   call solver(A, b, solution)

   deallocate (A, b)

   call info_print("PROGRAM END.")
end program main
