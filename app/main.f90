program main

   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE_2D

   implicit none

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
   !> basis_type: the type of the FE.
   integer :: basis_type = 202
      !! 201: 2D linear
      !! 202: 2D quadratic
   integer :: mesh_type = 504
      !! 503: 2D triangular mesh;
      !! 504: 2D quadrilateral mesh;

   real(wp), allocatable :: M(:, :)
   integer, allocatable :: T(:, :)
   integer, allocatable :: boundarynodes(:, :)
   integer, allocatable :: boundaryedges(:, :)

   ! real(wp), allocatable :: solution(:, :)

   ! integer :: N_basis
   ! real(wp) :: h_basis, max_FE_error

   ! integer :: index, index_2

   select case (basis_type)
   case (201)
      Nh_basis = Nh_parition
      Nv_basis = Nv_parition
   case (202)
      Nh_basis = 2*Nh_parition
      Nv_basis = 2*Nv_parition
   end select

   call info_print("PROGRAM START.")

   call field_info%init(left, right, bottom, top, &
                        Nh_parition, Nv_parition, &
                        Nh_basis, Nv_basis)

   call field_info%print('(A,2F6.3)')

   call generate_info_matrix(field_info, M, T, verbose=.true.)
   ! print *, "The information matrix generated."

   call generate_boundarynodes(field_info, boundarynodes, verbose=.true.)
   call generate_boundaryedges(field_info, boundaryedges, verbose=.true.)

   ! print *, "The boundary nodes and edges generated."

   call info_print("PROGRAM END.")
end program main
