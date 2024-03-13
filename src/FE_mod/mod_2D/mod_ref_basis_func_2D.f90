module mod_ref_basis_func_2D
   use mod_kinds
   implicit none
   private

   public :: triangular_ref_basis_func_2D

contains
   pure function triangular_ref_basis_func_2D(x, y, basis_type, basis_index, der_x, der_y) result(func)
      real(wp), intent(in) :: x, y
      integer, intent(in) :: basis_type, basis_index, der_x, der_y
      real(wp) :: func

      select case (basis_type)

      case (201) ! Linear basis function
         select case (der_x)
         case (0)
            select case (der_y)
            case (0)
               select case (basis_index)
               case (1)
                  func = 1.0_wp - x - y
               case (2)
                  func = x
               case (3)
                  func = y
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for linear basis function"
               end select
            case (1)
               select case (basis_index)
               case (1)
                  func = -1.0_wp
               case (2)
                  func = 0.0_wp
               case (3)
                  func = 1.0_wp
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for linear basis function"
               end select
            end select
         case (1)
            select case (der_y)
            case (0)
               select case (basis_index)
               case (1)
                  func = -1.0_wp
               case (2)
                  func = 1.0_wp
               case (3)
                  func = 0.0_wp
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for linear basis function"
               end select
            !! case (1) der_x = der_y = 1 is not implemented.
            case default
               error stop "triangular_ref_basis_func_2D: Invalid basis index for linear basis function, der_x = 1, der_y /= 0"
            end select
         end select

      case (202) ! Quadratic basis function
         select case (der_x)
         case (0)
            select case (der_y)
            case (0)
               select case (basis_index)
               case (1)
                  func = 1.0_wp - 3.0_wp*(x + y) + 2.0_wp*(x + y)**2
               case (2)
                  func = 2.0_wp*x**2 - x
               case (3)
                  func = 2.0_wp*y**2 - y
               case (4)
                  func = 4.0_wp*x*(1.0_wp - x - y)
               case (5)
                  func = 4.0_wp*x*y
               case (6)
                  func = 4.0_wp*y*(1.0_wp - x - y)
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for quadratic basis function"
               end select
            case (1)
               select case (basis_index)
               case (1)
                  func = -3.0_wp + 4.0_wp*(x + y)
               case (2)
                  func = 0.0_wp
               case (3)
                  func = 4.0_wp*y - 1.0_wp
               case (4)
                  func = -4.0_wp*x
               case (5)
                  func = 4.0_wp*x
               case (6)
                  func = 4.0_wp - 8.0_wp*y - 4.0_wp*x
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for quadratic basis function"
               end select
            case (2)
               select case (basis_index)
               case (1)
                  func = 4.0_wp
               case (2)
                  func = 0.0_wp
               case (3)
                  func = 4.0_wp
               case (4)
                  func = 0.0_wp
               case (5)
                  func = 0.0_wp
               case (6)
                  func = -8.0_wp
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for quadratic basis function"
               end select
            case default
               error stop "triangular_ref_basis_func_2D: Invalid der_y for quadratic basis function"
            end select
         case (1)
            select case (der_y)
            case (0)
               select case (basis_index)
               case (1)
                  func = -3.0_wp + 4.0_wp*(x + y)
               case (2)
                  func = 4.0_wp*x - 1.0_wp
               case (3)
                  func = 0
               case (4)
                  func = 4.0_wp - 8.0_wp*x - 4.0_wp*y
               case (5)
                  func = 4.0_wp*y
               case (6)
                  func = -4.0_wp*y
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for quadratic basis function"
               end select
            case (1)
               select case (basis_index)
               case (1)
                  func = 4.0_wp
               case (2)
                  func = 0.0_wp
               case (3)
                  func = 0.0_wp
               case (4)
                  func = -4.0_wp
               case (5)
                  func = 4.0_wp
               case (6)
                  func = -4.0_wp
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for quadratic basis function"
               end select
            case default
               error stop "triangular_ref_basis_func_2D: Invalid der_y for quadratic basis function"
            end select
         case (2)
            select case (der_y)
            case (0)
               select case (basis_index)
               case (1)
                  func = 4.0_wp
               case (2)
                  func = 4.0_wp
               case (3)
                  func = 0.0_wp
               case (4)
                  func = -8.0_wp
               case (5)
                  func = 0.0_wp
               case (6)
                  func = 0.0_wp
               case default
                  error stop "triangular_ref_basis_func_2D: Invalid basis index for quadratic basis function"
               end select
            case default
               error stop "triangular_ref_basis_func_2D: Invalid der_y for quadratic basis function"
            end select
         case default
            error stop "triangular_ref_basis_func_2D: Invalid der_x for quadratic basis function"
         end select

      end select

   end function triangular_ref_basis_func_2D

end module mod_ref_basis_func_2D
