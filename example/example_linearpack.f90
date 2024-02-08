program example_linearpack
! use logger_mod, only: logger_init, logger => master_logger
   use square_matrix_mod, only: dp
   use general_square_matrix_mod, only: general_square_matrix

   implicit none

   type(general_square_matrix) :: solver
   integer, parameter :: n = 5
   real(dp), dimension(n, n) :: matrix
   real(dp), dimension(n) :: x_actual, x_solved, b
   integer :: i

!> flogging
!>
! Initialise the logger prior to use
! call logger_init('log.out')

! Write some debugging information
! call logger%debug('logger_example', 'Starting program logger_example')

! Perform some calculation
! ...
! call logger%info('logger_example', 'Found result of calculation')

! Perform another calculation
! ...
! Oh no, an error has occurred
! call logger%error('logger_example', 'Calculation failed due to error')

! call logger%debug('logger_example', 'Ending program logger_example')

!> lapack-wrapper
!>
   matrix(1, 1) = 3.5_dp
   matrix(1, 2) = 1._dp
   matrix(1, 3) = -5._dp
   matrix(1, 4) = 1._dp
   matrix(1, 5) = 0._dp

   matrix(2, 1) = -0.5_dp
   matrix(2, 2) = 0.03_dp
   matrix(2, 3) = 8._dp
   matrix(2, 4) = 0._dp
   matrix(2, 5) = -7._dp

   matrix(3, 1) = -2.2_dp
   matrix(3, 2) = 100._dp
   matrix(3, 3) = 0._dp
   matrix(3, 4) = -1._dp
   matrix(3, 5) = -1._dp

   matrix(4, 1) = 5.5_dp
   matrix(4, 2) = 0._dp
   matrix(4, 3) = -11._dp
   matrix(4, 4) = -82._dp
   matrix(4, 5) = 2._dp

   matrix(5, 1) = 0._dp
   matrix(5, 2) = 5._dp
   matrix(5, 3) = 4._dp
   matrix(5, 4) = 3._dp
   matrix(5, 5) = -6._dp

   write (*, "(A)") "Solving linear system"
   do i = 1, n
      write (*, "(5F9.2)") matrix(i, :)
   end do

   solver = general_square_matrix(matrix)

   x_actual(1) = 1._dp
   x_actual(2) = 2._dp
   x_actual(3) = 3._dp
   x_actual(4) = 4._dp
   x_actual(5) = 5._dp

   b(1) = -5.5_dp
   b(2) = -11.44_dp
   b(3) = 188.8_dp
   b(4) = -345.5_dp
   b(5) = 7._dp

   write (*, "(/, A, T25, '[', 5F9.2, ']')") "RHS of linear system is", b
   write (*, "(A, T25, '[', 5F9.2, ']')") "Expected solution is", x_actual

   b = solver*x_actual
   x_solved = solver%inv_mat_mult(b)

   write (*, "(/, A, T25, '[', 5F9.2, ']')") "Actual solution is", x_solved
   write (*, "(A, F8.3)") "Backward Error = ", sqrt(sum((solver*x_solved - b)**2))
   write (*, "(A, F8.3)") "Forward Error  = ", sqrt(sum((x_solved - x_actual)**2))

end program example_linearpack
