module general_square_matrix_mod
   use square_matrix_mod
   use quadrature_module, only: wp => quadrature_wp
   implicit none
   private

   real(wp), dimension(:), allocatable :: real_work_array
   integer, dimension(:), allocatable :: integer_work_array

   type, extends(square_matrix), public :: general_square_matrix
      private
      real(wp), dimension(:, :), allocatable :: matrix
      real(wp), dimension(:, :), allocatable :: factored_matrix
      real(wp), dimension(:), allocatable :: col_scale_factors
      real(wp), dimension(:), allocatable :: row_scale_factors
      integer, dimension(:), allocatable :: pivots
      logical :: factored = .false.
      character :: equilibration = "N"
   contains
      procedure :: get_matrix
      procedure :: get_size
      procedure :: inv_mat_mult
   end type general_square_matrix

   interface general_square_matrix
      module procedure constructor
   end interface general_square_matrix

   interface
      pure subroutine dgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, &
                             ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, &
                             iwork, info)
         import wp
         character, intent(in) :: fact
         character, intent(in) :: trans
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         real(wp), intent(inout), dimension(lda, n) :: a
         integer, intent(in) :: lda
         real(wp), intent(inout), dimension(ldaf, n) :: af
         integer, intent(in) :: ldaf
         integer, intent(inout), dimension(n) :: ipiv
         character, intent(inout) :: equed
         real(wp), intent(inout), dimension(n) :: r
         real(wp), intent(inout), dimension(n) :: c
         real(wp), intent(inout), dimension(ldb, nrhs) :: b
         integer, intent(in) :: ldb
         real(wp), intent(out), dimension(ldx, nrhs) :: x
         integer, intent(in) :: ldx
         real(wp), intent(out) :: rcond
         real(wp), intent(out), dimension(nrhs) :: ferr
         real(wp), intent(out), dimension(nrhs) :: berr
         real(wp), intent(out), dimension(4*n) :: work
         integer, intent(out), dimension(n) :: iwork
         integer, intent(out) :: info
      end subroutine dgesvx
   end interface

contains

   function constructor(matrix) result(this)
      real(wp), dimension(:, :), intent(in) :: matrix
      type(general_square_matrix) :: this
      integer :: n
      n = size(matrix, 1)
      if (n /= size(matrix, 2)) then
         error stop "Non-square matrix used to build type(general_square_matrix)"
      end if
      allocate (this%matrix(n, n), source=matrix)
      allocate (this%factored_matrix(n, n))
      allocate (this%col_scale_factors(n))
      allocate (this%row_scale_factors(n))
      allocate (this%pivots(n))
      this%col_scale_factors = 1.0_dp
      this%row_scale_factors = 1.0_dp
   end function constructor

   pure function get_matrix(this) result(matrix)
      class(general_square_matrix), intent(in) :: this
      real(wp), dimension(this%get_size(), this%get_size()) :: matrix
      integer :: i, j
      do concurrent(i=1:this%get_size(), j=1:this%get_size())
         matrix(i, j) = this%matrix(i, j)
         if (this%equilibration == 'R' .or. this%equilibration == 'B') &
            matrix(i, j) = matrix(i, j)/this%row_scale_factors(i)
         if (this%equilibration == 'C' .or. this%equilibration == 'B') &
            matrix(i, j) = matrix(i, j)/this%col_scale_factors(j)
      end do
   end function get_matrix

   pure function get_size(this) result(mat_size)
      class(general_square_matrix), intent(in) :: this
      integer :: mat_size
      mat_size = size(this%matrix, 1)
   end function get_size

   function inv_mat_mult(this, rhs) result(solution)
      class(general_square_matrix), intent(inout) :: this
      real(wp), dimension(:), intent(in) :: rhs
      real(wp), dimension(size(rhs)) :: solution

      integer :: n, info
      real(wp), dimension(size(rhs), 1) :: x, b
      real(wp), dimension(1) :: ferr, berr
      real(wp) :: rcond
      character :: fact

      n = size(rhs)
      if (n /= this%get_size()) error stop "Mismatched vector/matrix sizes"

      if (.not. allocated(real_work_array)) then
         allocate (real_work_array(4*n))
         allocate (integer_work_array(n))
      else
         if (size(real_work_array) < 4*n) then
            deallocate (real_work_array)
            allocate (real_work_array(4*n))
         end if
         if (size(integer_work_array) < n) then
            deallocate (integer_work_array)
            allocate (integer_work_array(n))
         end if
      end if
      b(:, 1) = rhs

      if (this%factored) then
         fact = 'F'
      else
         fact = 'E'
      end if
      call dgesvx(fact, 'N', n, 1, this%matrix, n, this%factored_matrix, &
                  n, this%pivots, this%equilibration, this%row_scale_factors, &
                  this%col_scale_factors, b, n, x, n, rcond, ferr, berr, &
                  real_work_array, integer_work_array, info)

      if (info /= 0) error stop "DGESVX returned with non-zero INFO argument"
      this%factored = .true.
      solution = x(:, 1)
   end function inv_mat_mult

end module general_square_matrix_mod
