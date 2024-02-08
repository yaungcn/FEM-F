program example_FE_2D
   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE_2D

   implicit none

   !> input information of field
   type(field) :: input_info

   !> Field parameter input.
   real(wp) :: left = 0.0_wp, right = 1.0_wp
   !! The problem domain is [left,right]*[bottom,top].
   real(wp) :: h_partition = 1.0_wp/4.0_wp
   !! The step size of the partition.
   integer :: Gauss_point_number = 8
   !! The number of Gauss Quadrature points.
   !> basis_type: the type of the FE.
   integer :: basis_type = 101
   !! 101: 1D linear
   !! 102: undefined

   real(wp), allocatable :: solution(:, :)

   integer :: N_basis
   real(wp) :: h_basis, max_FE_error

   integer :: index

   print *, "------------------start------------------"

   do index = 1, 8
      !> Initial input information
      call input_info%init(left, right, h_partition, Gauss_point_number, basis_type)

      call fe_solver(solution, input_info)

      ! print *, solution

      if (input_info%basis_type == 101) then
         N_basis = int((input_info%right - input_info%left)/h_partition)
         h_basis = h_partition
      end if
      max_FE_error = max_FE_err(solution, N_basis, input_info%left, h_basis)
      write (*, '(A,E15.5)') '   h_partition = ', h_partition
      write (*, '(A,E15.5)') '    Max error  = ', max_FE_error
      print *, "-----------------------------------------"

      deallocate (solution)
      h_partition = 1.0_wp/(2.0_wp)**(index + 2)
   end do

end program example_FE_2D

