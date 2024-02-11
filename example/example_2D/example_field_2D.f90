program example_field_2D
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

   select case (basis_type)
   case (201)
      Nh_basis = Nh_parition
      Nv_basis = Nv_parition
   case (202)
      Nh_basis = 2*Nh_parition
      Nv_basis = 2*Nv_parition
   end select

   print *, "------------------Field Information------------------"
   !! TODO: Input the field information.
   !> Initialize and check the field information.
   call field_info%init(left, right, bottom, top, &
                        Nh_parition, Nv_parition, &
                        Nh_basis, Nv_basis, &
                        Gauss_point_number, &
                        tolerance=1.0e-6_wp, &
                        basis_type=basis_type, &
                        mesh_type=mesh_type)

   !> Use the format (A,2F6.3) to print the field information.
   call field_info%print('(A,2F6.3)')

   print *, "------------------Field Information-----------------"
end program example_field_2D
