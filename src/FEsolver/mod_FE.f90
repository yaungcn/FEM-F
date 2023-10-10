module mod_FE

   use mod_error
   use mod_field
   use mod_co_func
   use mod_assemble
   use mod_generate
   use mod_gauss_quad
   use mod_basis_func
   use mod_treat_boundary
   use mod_linear_solver
   use quadrature_module, wp => quadrature_wp

contains

   subroutine fe_solver(solution, field_info)
      real(wp), dimension(:, :), allocatable :: solution
      !> Parameter in FE solver
      type(field), intent(in) :: field_info

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
      ! real(wp), allocatable :: solution(:, :)

      type(func_a) :: cofunc_a
      type(func_f) :: cofunc_f
      type(func_g) :: cofunc_g

      integer :: trial_basis_type = 101, trial_derivate_degree = 1
      integer :: test_basis_type = 101, test_derivate_degree = 1, b_test_derivate_degree = 0

      call generate_info_matrix(M=M_partition, T=T_partition, field_info=field_info)
      !! Mesh information for partition and FE basis functions.

      N_partition = int((field_info%right - field_info%left)/field_info%h_partition)
      num_of_elements = N_partition

      select case (field_info%basis_type)
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
                              field_info)

      call assemble_vector_1D(b, cofunc_f, &
                              M_partition, T_partition, &
                              T_basis, num_of_elements, &
                              num_of_test_local_basis, &
                              test_basis_type, b_test_derivate_degree, &
                              field_info)

      boundarynodes = generate_boundarynodes(N_basis)
      !! Generate boundary nodes
      call treat_Dirchlet_boundary(cofunc_g, A, b, boundarynodes, M_basis)

      call solver(A, b, solution)
      deallocate (A, b)
   end subroutine
end module mod_FE
