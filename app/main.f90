program main
   use logger_mod, only: logger_init, logger => master_logger
   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE

   implicit none

   !> input information of field
   type(field) :: input_info

   !> Field parameter input.
   real(wp) :: left = 0.0_wp, right = 1.0_wp
         !! The problem domain is [left,right]*[bottom,top].
   real(wp) :: h_partition = 1.0_wp/256.0_wp
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

   !> flogging
   call logger_init('./log/log.out')
      !! Initialise the logger prior to use
   call logger%info('main_log', 'Program Starts')
      !! log information

   !> Initial input information
   call input_info%init(left, right, h_partition, Gauss_point_number, basis_type)

   call fe_solver(solution, input_info)

   ! print *, solution

   if (input_info%basis_type == 101) then
      N_basis = int((input_info%right - input_info%left)/h_partition)
      h_basis = h_partition
   end if
   max_FE_error = max_FE_err(solution, N_basis, input_info%left, h_basis)
   print *, 'h_partition = ', h_partition
   print *, 'Max error = ', max_FE_error

   call logger%info('main_log', 'Program Ends')

   deallocate (solution)
end program main

