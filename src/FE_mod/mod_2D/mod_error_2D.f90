module mod_error_2D
   use mod_kinds
   implicit none

contains
   pure function max_FE_err(solution, N_basis, left, h_basis)
      real(wp), intent(in) :: solution(:, :)
      integer, intent(in) :: N_basis
      real(wp), intent(in) :: left, h_basis
      real(wp) :: max_FE_err, err
      integer :: index, index2
      max_FE_err = 0
      do index = 1, N_basis + 1
         do index2 = 1, N_basis + 1
            err = solution(index, index2) - exact_solution(left + (index - 1)*h_basis, left + (index2 - 1)*h_basis)
            if (abs(max_FE_err) < abs(err)) then
               max_FE_err = err
            end if
         end do
      end do
   end function

   pure function exact_solution(x, y)
      !! f(x) = x*cos(x)
      real(wp), intent(in) :: x, y
      real(wp) :: exact_solution
      exact_solution = x*y*(1 - x/2)*(1 - y)*exp(x + y)
   end function
end module mod_error_2D
