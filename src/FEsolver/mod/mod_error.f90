module mod_error
   use quadrature_module, only: wp => quadrature_wp
   implicit none

contains
   function max_FE_err(solution, N_basis, left, h_basis)
      real(wp), intent(in) :: solution(:, :)
      integer, intent(in) :: N_basis
      real(wp), intent(in) :: left, h_basis
      real(wp) :: max_FE_err, err
      integer :: index
      max_FE_err = 0
      do index = 1, N_basis + 1
         err = solution(index, 1) - exact_solution(left + (index - 1)*h_basis)
         if (abs(max_FE_err) < abs(err)) then
            max_FE_err = err
         end if
      end do
   end function

   function exact_solution(x)
      !! f(x) = x*cos(x)
      real(wp), intent(in) :: x
      real(wp) :: exact_solution
      exact_solution = x*cos(x)
   end function
end module mod_error
