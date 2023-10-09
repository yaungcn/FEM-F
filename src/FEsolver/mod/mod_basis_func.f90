module mod_basis_func
   use quadrature_module, only: wp => quadrature_wp
   implicit none

   type, abstract :: local_basis
   contains
      procedure(basis_func), nopass, deferred :: basis_f
   end type local_basis

   abstract interface
      function basis_func(x, vertices, basis_type, basis_index, derivative_degree) &
         result(f)
         import :: wp
      !! This is for the local basis functions of 1D FE.
         real(wp), intent(in) :: x
            !! x: the coordinate of the point where we want to
            !!    evaluate the local FE basis function.
         real(wp), intent(in) :: vertices(2)
            !! vertices: start[vertice(1)] and end[vertice(2)] of the local partition
         integer, intent(in) :: basis_type, basis_index, derivative_degree
            !! basis_type: the type of the FE.
            !! basis_type=101:1D linear FE.
            !! basis_index: the index of basis function to specify
            !!              which basis function we want to use.
            !! derivative_degree:the derivative degree of the FE basis function.
         real(wp) :: f
      end function basis_func
   end interface

   type, extends(local_basis) :: local_basis_1D
   contains
      procedure, nopass :: basis_f => basis_func1D
   end type local_basis_1D

contains
   function basis_func1D(x, vertices, basis_type, basis_index, derivative_degree) &
      result(f)
      !! This is for the local basis functions of 1D FE.
      real(wp), intent(in) :: x
      real(wp), intent(in) :: vertices(2)
      integer, intent(in) :: basis_type, basis_index, derivative_degree
      real(wp) :: f

      select case (basis_type)
      case (101)

         select case (derivative_degree)
         case (0)
            select case (basis_index)
            case (1)
               f = (vertices(2) - x)/(vertices(2) - vertices(1))
            case (2)
               f = (x - vertices(1))/(vertices(2) - vertices(1))
            end select
         case (1)
            select case (basis_index)
            case (1)
               f = 1.0_wp/(vertices(1) - vertices(2))
            case (2)
               f = 1.0_wp/(vertices(2) - vertices(1))
            end select
         end select

      case default
         error stop "Invalid basis_type! "
      end select

   end function basis_func1D
end module mod_basis_func
