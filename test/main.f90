program main

   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use M_attr
   use mod_FE_2D

   implicit none
   abstract interface
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
   real(wp) :: left = 0.0_wp, right = 1.0_wp, bottom = 0.0_wp, top = 1.0_wp
      !! The problem domain is [left,right]*[bottom,top].
   integer :: Nh_parition = 2, Nv_parition = 2
      !! The number of partition in horizontal and vertical direction.
   ! real(wp) :: h_partition, v_partition
      !! The step size of the partition.
   integer :: Nh_basis, Nv_basis
      !! The number of FE basis functions in horizontal and vertical direction.
   integer :: Gauss_point_number = 8
      !! The number of Gauss Quadrature points.
   real(wp) :: tolerance = 1.0e-6_wp
      !! The tolerance of the error.
   !> basis_type: the type of the FE.
   integer :: basis_type = 201
      !! 201: 2D linear
      !! 202: 2D quadratic
   integer :: mesh_type = 503
      !! 503: 2D triangular mesh;
      !! 504: 2D quadrilateral mesh;

   real(wp), allocatable :: M(:, :), M_basis(:, :)
   integer, allocatable :: T(:, :), T_basis(:, :)
   integer, allocatable :: boundarynodes(:, :)
   integer, allocatable :: boundaryedges(:, :)

   character(len=4096)   :: arg0

   real(wp), allocatable :: A(:, :), b(:, :), A1(:, :), A2(:, :)
   real(wp), allocatable :: solution(:, :)
   real(wp) :: vet(2, 3)

   procedure(local_basis_func), pointer :: lbfunc => null()
   procedure(cofunc_2d), pointer :: cofunc => null()
   procedure(cofunc_2d), pointer :: cofuncf => null()
   lbfunc => local_basis_func_2D
   cofunc => function_c
   cofuncf => function_f

   select case (basis_type)
   case (201)
      Nh_basis = Nh_parition
      Nv_basis = Nv_parition
   case (202)
      Nh_basis = 2*Nh_parition
      Nv_basis = 2*Nv_parition
   end select

   call attr_mode(manner='color') !! plain, color, raw, ansi
   call info_print("PROGRAM START.")

   call get_command_argument(0, arg0)
   if (index(arg0, '/') .ne. 0) arg0 = arg0(index(arg0, '/', back=.true.) + 1:)
   if (index(arg0, '\') .ne. 0) arg0 = arg0(index(arg0, '\', back=.true.) + 1:)
   call alert("DATE   ", "<YE><MO><DA> <HR>:<MI>:<SE>.<MS> <TZ>", arg0)

   call field_info%init(left, right, bottom, top, &
                        Nh_parition, Nv_parition, &
                        Nh_basis, Nv_basis, &
                        Gauss_point_number=Gauss_point_number, &
                        tolerance=tolerance, &
                        trial_basis_type=basis_type, &
                        test_basis_type=basis_type, &
                        mesh_type=mesh_type)

   call field_info%print('(A,2F6.3)')

   call generate_info_matrix(field_info, M, T, verbose=.true.)
   !!!
   print *, "M = ", M
   print *, "T = ", T
   vet(:, :) = M(:, T(:, 1))
   print *, "vet = ", vet
   print *, lbfunc(0.5_wp, 0.5_wp, vet, 1, 1, 0, 0)
   !!!
   call generate_boundarynodes(field_info, boundarynodes, verbose=.true.)
   call generate_boundaryedges(field_info, boundaryedges, verbose=.true.)

   select case (basis_type)
   case (201)
      M_basis = M
      T_basis = T
   case (202)
      error stop "Not implemented yet."
   end select

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

   call assemble_vector_2D(b, cofuncf, lbfunc, &
                           M, T, &
                           T_basis, &
                           1, 1, &
                           field_info)

   call treat_Dirchlet_boundary(function_g, A, b, boundarynodes, M)

   call solver(A, b, solution)

   call info_print("PROGRAM END.")
end program main
