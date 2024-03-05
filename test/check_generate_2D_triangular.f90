program check_generate_2D_triangular
   use quadrature_module, wp => quadrature_wp
   use mod_field_2D
   use mod_generate_2D

   implicit none

   !> input information of field
   type(field) :: field_info

   !> Field parameter input.
   real(wp) :: left = 0.0_wp, right = 1.0_wp, bottom = 0.0_wp, top = 1.0_wp
      !! The problem domain is [left,right]*[bottom,top].
   integer :: Nh = 8, Nv = 8
      !! The number of partition in horizontal and vertical direction.
   integer :: Nh_basis, Nv_basis
      !! The number of FE basis functions in horizontal and vertical direction.
   integer :: Gauss_point_number = 8
      !! The number of Gauss Quadrature points.
   integer :: basis_type = 201
      !! 201: 2D linear
      !! 202: 2D quadratic
   integer :: mesh_type = 503
      !! 503: 2D triangular mesh;
      !! 504: 2D quadrilateral mesh;

   !> infortramtion matrix, M and T
   real(wp), allocatable :: M(:, :)
   integer, allocatable :: T(:, :)
   integer :: index

   select case (basis_type)
   case (201)
      Nh_basis = Nh
      Nv_basis = Nv
   case (202)
      Nh_basis = 2*Nh
      Nv_basis = 2*Nv
   end select

   call field_info%init(left, right, bottom, top, &
                        Nh, Nv, &
                        Nh_basis, Nv_basis, &
                        Gauss_point_number=Gauss_point_number, &
                        tolerance=1.0e-6_wp, &
                        trial_basis_type=basis_type, &
                        test_basis_type=basis_type, &
                        mesh_type=mesh_type)
   call field_info%check()
   write (*, '(A,I3)') "mesh_type = ", field_info%mesh_type

   call generate_info_matrix(field_info, M, T)

   write (*, '(A,I3.3,A,I3.3)') "the mesh size: ", Nv, "x", Nh
   write (*, '(A,I1.1,A,I3.3)') "M = ", size(M, 1), "x", size(M, 2)
   do index = 1, size(M, 1)
      print *, M(index, :)
   end do

   write (*, '(A,I1.1,A,I3.3)') "T = ", size(T, 1), "x", size(T, 2)
   do index = 1, size(T, 1)
      print *, T(index, :)
   end do

   print *, "noe=", field_info%number_of_elements

end program check_generate_2D_triangular
