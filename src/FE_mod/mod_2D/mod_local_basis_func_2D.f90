module mod_local_basis_func_2D
   use quadrature_module, only: wp => quadrature_wp
   use mod_ref_basis_func_2D, only: trf => triangular_ref_basis_func_2D
   implicit none
   private

   public :: local_basis_func_2D

contains
   pure function local_basis_func_2D(x, y, vertices, basis_type, basis_index, der_x, der_y) result(func)
      !! the local basis function is generated from the reference basis function
      real(wp), intent(in) :: x, y
      real(wp), intent(in) :: vertices(2, 3)
      real(wp) :: x_hat, y_hat, J_11, J_12, J_21, J_22, J_det
      integer, intent(in) :: basis_type, basis_index, der_x, der_y
      real(wp) :: func

      J_11 = vertices(1, 2) - vertices(1, 1)
      J_12 = vertices(1, 3) - vertices(1, 1)
      J_21 = vertices(2, 2) - vertices(2, 1)
      J_22 = vertices(2, 3) - vertices(2, 1)
      J_det = J_11*J_22 - J_12*J_21

      x_hat = (J_22*(x - vertices(1, 1)) - J_12*(y - vertices(2, 1)))/J_det
      y_hat = (-J_21*(x - vertices(1, 1)) + J_11*(y - vertices(2, 1)))/J_det
      select case (der_x)
      case (0)
         select case (der_y)
         case (0)
            func = trf(x_hat, y_hat, basis_type, basis_index, 0, 0)
         case (1)
            func = (trf(x_hat, y_hat, basis_type, basis_index, 1, 0)*(-J_12) + &
                    trf(x_hat, y_hat, basis_type, basis_index, 0, 1)*J_11)/J_det
         case (2)
            func = (trf(x_hat,y_hat,basis_type,basis_index,2,0)*J_12**2+trf(x_hat,y_hat,basis_type,basis_index,0,2)*J_11**2+trf(x_hat,y_hat,basis_type,basis_index,1,1)*(-2*J_11*J_12))/J_det**2
         case default
            error stop 'Error: der_y is not valid'
         end select
      case (1)
         select case (der_y)
         case (0)
            func = (trf(x_hat, y_hat, basis_type, basis_index, 1, 0)*J_22 + &
                    trf(x_hat, y_hat, basis_type, basis_index, 0, 1)*(-J_21))/J_det
         case (1)
            func = (trf(x_hat, y_hat, basis_type, basis_index, 2, 0)*(-J_22*J_12) + &
                    trf(x_hat, y_hat, basis_type, basis_index, 0, 2)*(-J_21*J_11) + &
                    trf(x_hat, y_hat, basis_type, basis_index, 1, 1)*(J_21*J_12 + J_11*J_22))/J_det**2
         case default
            error stop 'Error: der_y is not valid'
         end select
      case (2)
         select case (der_y)
         case (0)
            func = (trf(x_hat, y_hat, basis_type, basis_index, 2, 0)*J_12**2 + &
                    trf(x_hat, y_hat, basis_type, basis_index, 0, 2)*J_11**2 + &
                    trf(x_hat, y_hat, basis_type, basis_index, 1, 1)*(-2*J_11*J_12)) &
                   /J_det**2
         case default
            error stop 'Error: der_y is not valid'
         end select
      end select

   end function local_basis_func_2D
end module mod_local_basis_func_2D
