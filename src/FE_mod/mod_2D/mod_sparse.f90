module mod_sparse
   use mod_kinds
   use lsqr_module, only: lsqr_solver_ez
   implicit none
   private

   public :: sparse_COO

   !> The sparse matrix type
   type sparse_COO
      !> The number of columns and rows in the matrix
      integer :: nrow, ncol
      integer :: num_nonzeros = 0
      !> The column and row index of the non-zero values
      integer, dimension(:), allocatable :: icol, irow
      !> The non-zero values
      real(wp),dimension(:), allocatable :: values

   contains
      procedure, pass(self) :: init => init_sparse
      procedure, pass(self) :: append => append_nonzero
      procedure, pass(self) :: to_dense
      procedure, pass(self) :: from_dense
      procedure, pass(self) :: solver => solver_sparse
      procedure, pass(self) :: destory => finalize_sparse
      procedure, pass(self) :: finalize => finalize_sparse

      procedure, pass(self) :: plus_sparse
      generic :: operator(+) => plus_sparse  !! overload the `+` operator
   end type sparse_COO
contains
   !> Initialize the sparse matrix space
   pure subroutine init_sparse(self, nrow, ncol)

      class(sparse_COO), intent(inout) :: self
      integer, intent(in) :: nrow, ncol

      self%nrow=nrow
      self%ncol=ncol

      allocate(self%icol(self%num_nonzeros), self%irow(self%num_nonzeros), self%values(self%num_nonzeros))

   end subroutine init_sparse
   !> Append a value to the sparse matrix
   pure subroutine append_nonzero(self, row_index, col_index, value_append)

      class(sparse_COO), intent(inout) :: self

      !> The row and column index of the value to be appended
      integer, intent(in) :: row_index, col_index
      !> The value to be appended
      real(wp), intent(in) :: value_append

      ! Increase the number of non-zero elements
      self%num_nonzeros = self%num_nonzeros + 1

      ! Resize the arrays
      call resize_index(self%irow, self%num_nonzeros)
      call resize_index(self%icol, self%num_nonzeros)
      call resize_value(self%values, self%num_nonzeros)

      ! Add the new non-zero element
      self%irow(self%num_nonzeros) = row_index
      self%icol(self%num_nonzeros) = col_index
      self%values(self%num_nonzeros) = value_append

   end subroutine append_nonzero
   !> Resize an integer array
   pure subroutine resize_index(array, new_size)
      implicit none
      integer, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: new_size
      integer, allocatable :: temp(:)

      ! Allocate a temporary array with the new size
      allocate(temp(new_size))

      ! Copy the old elements to the temporary array
      if (allocated(array)) then
         temp(1:size(array)) = array
         deallocate(array)
      end if

      ! Rename the temporary array
      call move_alloc(temp, array)
   end subroutine resize_index
   !> Resize a real array
   pure subroutine resize_value(array, new_size)
      implicit none
      real(wp), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: new_size
      real(wp), allocatable :: temp(:)

      ! Allocate a temporary array with the new size
      allocate(temp(new_size))

      ! Copy the old elements to the temporary array
      if (allocated(array)) then
         temp(1:size(array)) = array
         deallocate(array)
      end if

      ! Rename the temporary array
      call move_alloc(temp, array)
   end subroutine resize_value
   !> Convert the sparse matrix to a dense matrix
   pure subroutine to_dense(self, dense_mat)
      class(sparse_COO), intent(in) :: self
      real(wp), dimension(:,:), allocatable, intent(out) :: dense_mat
      integer :: i

      ! initialize the dense matrix
      allocate(dense_mat(self%nrow, self%ncol))
      dense_mat = 0.0_wp

      ! copy the non-zero values to the dense matrix
      do i = 1, self%num_nonzeros
         dense_mat(self%irow(i), self%icol(i)) = self%values(i)
      end do

   end subroutine to_dense
   !> Convert the dense matrix to a sparse matrix
   pure subroutine from_dense(self, dense_mat)
      class(sparse_COO), intent(inout) :: self
      real(wp), dimension(:,:), intent(in) :: dense_mat
      integer :: i, j

      ! initialize the sparse matrix
      call self%init(size(dense_mat, 1), size(dense_mat, 2))

      ! loop over the dense matrix and append the non-zero values
      do i = 1, size(dense_mat, 1)
         do j = 1, size(dense_mat, 2)
            if (dense_mat(i, j) /= 0.0_wp) then
               call self%append(i, j, dense_mat(i, j))
            end if
         end do
      end do

   end subroutine from_dense
   !> Solver of the sparse matrix
   subroutine solver_sparse(self, b, x, istop)
      class(sparse_COO), intent(inout) :: self
      type(lsqr_solver_ez) :: solver  !! main solver class
      real(wp),dimension(:), intent(in) :: b
      !! right-hand side of the linear system, `A*x = b`.
      !! `b` is not a sparse vector
      real(wp),dimension(:), intent(out) :: x
      !! solution to `A*x = b`
      integer, intent(out) :: istop  !! solver exit code

      !> Solve the sparse matrix
      call solver%initialize(self%nrow,self%ncol,self%values,self%irow,self%icol)
      !! use defaults for other optional inputs
      call solver%solve(b,zero,x,istop)
      !! solve the linear system
   end subroutine solver_sparse
   !> Finalize-destory the sparse matrix
   subroutine finalize_sparse(self)
      class(sparse_COO), intent(inout) :: self

      deallocate(self%icol, self%irow, self%values)

   end subroutine finalize_sparse
   !> Overload the `+` operator
   function plus_sparse(self, val) result(ans)
      class(sparse_COO), intent(in) :: self
      class(sparse_COO), intent(in) :: val
      type(sparse_COO) :: ans
      integer :: i,j,idx

      if (self%nrow /= val%nrow .or. self%ncol /= val%ncol) then
         error stop "Error: The dimensions of the two matrices are not the same"
      end if

      ans = self

      do i = 1, val%num_nonzeros
         idx = -1
         do j = 1, self%num_nonzeros
            if (self%irow(j) == val%irow(i) .and. self%icol(j) == val%icol(i)) then
               idx = j
               exit
            end if
         end do
         if (idx > 0) then
            ans%values(idx) = ans%values(idx) + val%values(i)
         else
            call ans%append(val%irow(i), val%icol(i), val%values(i))
         end if
      end do

   end function   plus_sparse
end module mod_sparse
