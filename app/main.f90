program main

   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE_2D

   implicit none

   !> input information of field
   type(field) :: field_info

   !> Field parameter input.
   real(wp) :: left = 0.0_wp, right = 1.0_wp, bottom = 0.0_wp, top = 1.0_wp
   !! The problem domain is [left,right]*[bottom,top].
   integer :: Nh = 2, Nv = 2
   !! The number of partition in horizontal and vertical direction.
   ! real(wp) :: h_partition, v_partition
   !! The step size of the partition.
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

   ! real(wp), allocatable :: solution(:, :)

   ! integer :: N_basis
   ! real(wp) :: h_basis, max_FE_error

   integer :: index, index_2

   print *, "------------------start------------------"

   call field_info%init(left, right, bottom, top, Nh, Nv, Gauss_point_number, basis_type, mesh_type)
   call field_info%check()
   write (*, '(A,I3)') "mesh_type = ", field_info%mesh_type

   call generate_info_matrix(M, T, field_info)
100 format(A, I2, A, I2, A)
   write (*, 100) "the mesh size: [", Nv, "x", Nh, ']'
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

end program main
