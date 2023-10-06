program main
   use logger_mod, only: logger_init, logger => master_logger
   use linear_pack, only: general_square_matrix
   use quadrature_module, only: wp => quadrature_wp
   use mod_FE

   implicit none

   real(wp), allocatable :: P, T
   ! mesh node information matrix
   real(wp), allocatable :: Pb, Tb
   ! finite element node information matrix
   integer, allocatable :: boundarynodes
   real(wp), allocatable :: A, b
   real(wp), allocatable :: solution
   integer :: basis_type_trial = 101
   !> 101: 1D linear
   !> 102:

   !> flogging
   call logger_init('./log/log.out')
   ! Initialise the logger prior to use
   call logger%info('main_log', 'Program Starts')
   ! log information

   call generate_info_matrix(P, T)

   select case (basis_type_trial)
   case (101)
      call generate_info_matrix(Pb, Tb)
   case (102)
      continue
   case default
      continue
   end select

   call generate_boundarynodes(boundarynodes)

   call assemble_matirx_1D(A)

   call assemble_vector_1D(b)

   call treat_Dirchlet_boundary(A, b)

   call solver(A, b, solution)

   call logger%info('main_log', 'Program Ends')
end program main

