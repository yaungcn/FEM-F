program example_generate_2D
         use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE_2D

   implicit none

   !> input information of field
   type(field) :: field_info

   !> Field parameter input.
   real(wp) :: left = 0.0_wp, right = 1.0_wp, bottom = 0.0_wp, top = 1.0_wp
      !! The problem domain is [left,right]*[bottom,top].
   integer :: Nh_parition = 2, Nv_parition = 2
      !! The number of partition in horizontal and vertical direction.
      !! real(wp) :: h_partition, v_partition
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


   integer :: index, index_2

   select case (basis_type)
   case (201)
      Nh_basis = Nh_parition
      Nv_basis = Nv_parition
   case (202)
      Nh_basis = 2*Nh_parition
      Nv_basis = 2*Nv_parition
   end select

   print *, "------------------Information Matrix Generate------------------"

   call field_info%init(left, right, bottom, top, &
                        Nh_parition, Nv_parition, &
                        Nh_basis, Nv_basis, &
                        Gauss_point_number, &
                        basis_type, mesh_type)
   !> Check the field information.
   call field_info%check()
   !> Print the field information.
   call field_info%print()

   !> Generate the information matrix.
   call generate_info_matrix(M, T, field_info)

100 format(A, I2.2, A, I2.2, A)
   write (*, 100) "the mesh size: [", Nv_parition, "x", Nh_parition, ']'
   write (*, 100) "M = [", size(M, 1), "x", size(M, 2), ']'
   do index = 1, size(M, 1)
      do index_2 = 1, size(M, 2)
         write (*, '(F6.2)', advance='no') M(index, index_2)
      end do
      print *
   end do

   write (*, 100) "T = [", size(T, 1), "x", size(T, 2), ']'
   do index = 1, size(T, 1)
      print *, T(index, :)
   end do

   boundarynodes = generate_boundarynodes(field_info)
   boundaryedges = generate_boundaryedges(field_info)

   write (*, 100) "boundarynodes = [", size(boundarynodes, 1), "x", size(boundarynodes, 2), ']'
   do index = 1, size(boundarynodes, 1)
      print *, boundarynodes(index, :)
   end do

   write (*, 100) "boundaryedges = [", size(boundaryedges, 1), "x", size(boundaryedges, 2), ']'
   do index = 1, size(boundaryedges, 1)
      print *, boundaryedges(index, :)
   end do

   print *, "------------------Information Matrix Generate------------------"
end program example_generate_2D