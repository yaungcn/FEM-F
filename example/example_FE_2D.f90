program example_FE_2D
   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE_2D

   implicit none

   !> input information of field
   type(field) :: field_info

   !> Field parameter input.
   real(wp) :: left = 0.0_wp, right = 1.0_wp, bottom = 0.0_wp, top = 1.0_wp
   !! The problem domain is [left,right]*[bottom,top].
   integer :: Nh = 3, Nv = 1
   !! The number of partition in horizontal and vertical direction.
   ! real(wp) :: h_partition, v_partition
   !! The step size of the partition.
   integer :: Gauss_point_number = 8
   !! The number of Gauss Quadrature points.
   !> basis_type: the type of the FE.
   integer :: basis_type = 201
   !! 201: 2D linear
   !! 202: 2D quadratic
   integer :: mesh_type = 503
   !! 503: 2D triangular mesh;
   !! 504: 2D quadrilateral mesh;

   real(wp), allocatable :: M(:, :)
   integer, allocatable :: T(:, :)

   real(wp), allocatable :: solution(:, :)

   integer :: N_basis
   real(wp) :: h_basis, max_FE_error

   integer :: index

   print *, "------------------start------------------"

   call field_info%init(left, right, bottom, top, Nh, Nv, Gauss_point_number, basis_type, mesh_type)
   call field_info%check()
   print *, "mesh_type = ", field_info%mesh_type

   call generate_info_matrix(M, T, field_info)

   print *, "M = "
   do index = 1, size(M, 1)
      print *, M(index, :)
   end do
   ! print *, M
   print *, "T = ", size(T, 1)
   do index = 1, size(T, 1)
      print *, T(index, :)
   end do
   ! print *, T

   ! do index = 1, 8
   !    !> Initial input information
   !    call field_info%init(left, right, bottom, top, Nh, Nv, Gauss_point_number, basis_type)

   !    call fe_solver(solution, field_info)

   !    ! print *, solution

   !    if (field_info%basis_type == 101) then
   !       N_basis = int((field_info%right - field_info%left)/h_partition)
   !       h_basis = h_partition
   !    end if
   !    max_FE_error = max_FE_err(solution, N_basis, field_info%left, h_basis)
   !    write (*, '(A,E15.5)') '   h_partition = ', h_partition
   !    write (*, '(A,E15.5)') '    Max error  = ', max_FE_error
   !    print *, "-----------------------------------------"

   !    deallocate (solution)
   !    h_partition = 1.0_wp/(2.0_wp)**(index + 2)
   ! end do

end program example_FE_2D

