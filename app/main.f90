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
   real(wp) :: h_partition = 1.0_wp/4.0_wp
         !! The step size of the partition.
   integer :: Gauss_point_number = 8
         !! The number of Gauss Quadrature points.
   !> basis_type: the type of the FE.
   integer :: basis_type = 101
      !! 101: 1D linear
      !! 102:

   !> Parameter in FE solver

   !> INFORMATION MATRIX
   !> M stores the coordinates of all nodes for the type of FE specified by "basis_type".
   !> -M(i,j) is the ith coordinate of the jth node.
   !> -In 1D, M has only one row.
   !> T stores the global indices of the nodes od every element for the type of FE specified by "basis_type".
   !> -T(i,j) stores the global index of the ith node in jth element.
   real(wp), allocatable, target :: M_basis(:, :), M_partition(:, :)
      !! FE nodes information matrix
   integer, allocatable :: T_basis(:, :), T_partition(:, :)
      !! Mesh nodes information matrix

   !> FE & Mesh setting.
   integer :: N_basis, N_partition, num_of_elements
      !! N_basis: The N for the FE basis functions, not the partition.
      !! N_partition: The N for the partition, not the FE basis functions.
      !! N: The number of elements(sub-intervals).
   integer :: num_of_trial_local_basis, num_of_test_local_basis

   integer, allocatable :: boundarynodes(:, :)

   real(wp), allocatable :: A(:, :), b(:, :)
   real(wp), allocatable :: solution(:, :)

   type(func_a) :: cofunc_a
   type(func_f) :: cofunc_f
   type(func_g) :: cofunc_g

   integer :: trial_basis_type = 101, trial_derivate_degree = 1
   integer :: test_basis_type = 101, test_derivate_degree = 1, b_test_derivate_degree = 0

   integer :: index
   real(wp) :: ans

   !> flogging
   call logger_init('./log/log.out')
      !! Initialise the logger prior to use
   call logger%info('main_log', 'Program Starts')
      !! log information

   !> Initial input information
   call input_info%init(left, right, h_partition, Gauss_point_number, basis_type)

   call generate_info_matrix(M=M_partition, T=T_partition, field_info=input_info)
      !! Mesh information for partition and FE basis functions.

   N_partition = int((input_info%right - input_info%left)/h_partition)
   num_of_elements = N_partition
   select case (basis_type)
   case (101)
      N_basis = N_partition

      M_basis = M_partition
      T_basis = T_partition

      num_of_trial_local_basis = 2
      num_of_test_local_basis = 2
   case (102)
      continue
   case default
      continue
   end select

   !> Allocate memory to A, b and solution. Set to zero.
   allocate (A(N_basis + 1, N_basis + 1), b(N_basis + 1, 1), solution(N_basis + 1, 1))
   A = 0
   b = 0
   solution = 0

   call assemble_matirx_1D(A, cofunc_a, &
                           M_partition, T_partition, &
                           T_basis, T_basis, &
                           num_of_elements, &
                           num_of_trial_local_basis, num_of_test_local_basis, &
                           trial_basis_type, trial_derivate_degree, &
                           test_basis_type, test_derivate_degree, &
                           input_info)

   call assemble_vector_1D(b, cofunc_f, &
                           M_partition, T_partition, &
                           T_basis, num_of_elements, &
                           num_of_test_local_basis, &
                           test_basis_type, b_test_derivate_degree, &
                           input_info)

! 09 format('(5(f15.5),A,f15.5)')
   ! print *, 'Assemble A & b'
   ! do index = 1, N_basis + 1
   !    write (*, '(5(f15.5),A,f15.5)') A(index, :), '  |', b(index, 1)
   ! end do

   boundarynodes = generate_boundarynodes(N_basis)
      !! Generate boundary nodes
   call treat_Dirchlet_boundary(cofunc_g, A, b, boundarynodes, M_basis)

   print *, 'Treated A & b'
   do index = 1, N_basis + 1
      write (*, '(5(f15.5),A,f15.5)') A(index, :), '  |', b(index, 1)
   end do
   call solver(A, b, solution)

   print *, 'Solution:'
   do index = 1, N_basis + 1
      write (*, '(5(f15.5))') solution(index, 1)
   end do

   call logger%info('main_log', 'Program Ends')

   deallocate (A, b, solution)
end program main

