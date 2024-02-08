module mod_basis_func_1D
   use quadrature_module, only: wp => quadrature_wp
   implicit none

   private
   public :: local_basis, local_basis_1D

   type, abstract :: local_basis
      real(wp) :: vertices(2)
            !! vertices: start[vertice(1)] and end[vertice(2)] of the local partition
      integer :: basis_type, basis_index, derivative_degree
            !! basis_type: the type of the FE.
            !! basis_type=101:1D linear FE.
            !! basis_index: the index of basis function to specify
            !!              which basis function we want to use.
            !! derivative_degree:the derivative degree of the FE basis function.
   contains
      procedure(basis_initial), pass(self), deferred :: init
      procedure(basis_func), pass(self), deferred :: basis_f
   end type local_basis

   abstract interface
      subroutine basis_initial(self, vertices, basis_type, basis_index, derivative_degree)
         import :: wp, local_basis
         class(local_basis), intent(inout) :: self
         real(wp), intent(in) :: vertices(2)
         integer, intent(in) :: basis_type, basis_index, derivative_degree
      end subroutine basis_initial

      function basis_func(self, x) &
         result(f)
         import :: wp, local_basis
      !! This is for the local basis functions of 1D FE.
         class(local_basis), intent(in) :: self
         real(wp), intent(in) :: x
            !! x: the coordinate of the point where we want to
            !!    evaluate the local FE basis function.
         real(wp) :: f
      end function basis_func
   end interface

   type, extends(local_basis) :: local_basis_1D
   contains
      procedure :: init => basis_init1D
      procedure :: basis_f => basis_func1D
   end type local_basis_1D

contains

   subroutine basis_init1D(self, vertices, basis_type, basis_index, derivative_degree)
      class(local_basis_1D), intent(inout) :: self
      real(wp), intent(in) :: vertices(2)
      integer, intent(in) :: basis_type, basis_index, derivative_degree

      self%vertices = vertices
      self%basis_type = basis_type
      self%basis_index = basis_index
      self%derivative_degree = derivative_degree

   end subroutine

   function basis_func1D(self, x) &
      result(f)
      !! This is for the local basis functions of 1D FE.
      class(local_basis_1D), intent(in) :: self
      real(wp), intent(in) :: x
      real(wp) :: f

      select case (self%basis_type)
      case (101)

         select case (self%derivative_degree)
         case (0)
            select case (self%basis_index)
            case (1)
               f = (self%vertices(2) - x)/(self%vertices(2) - self%vertices(1))
            case (2)
               f = (x - self%vertices(1))/(self%vertices(2) - self%vertices(1))
            end select
         case (1)
            select case (self%basis_index)
            case (1)
               f = 1.0_wp/(self%vertices(1) - self%vertices(2))
            case (2)
               f = 1.0_wp/(self%vertices(2) - self%vertices(1))
            end select
         end select
      case (102)
         error stop "Local basis type '102' undefined! "
      case default
         error stop "Invalid basis_type! "
      end select

   end function basis_func1D

end module mod_basis_func_1D
