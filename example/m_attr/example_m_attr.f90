program mattr_test
   use M_attr, only: attr, attr_mode, attr_update, alert
   implicit none
   real :: x, y, z
   character(len=256) :: line
   character(len=*), parameter :: f = '( &
    &"   <bo><w><G> GREAT: </G></w>&
    &The new value <Y><b>",1x,a,1x,"</b></Y> is in range"&
    &)'
   character(len=*), parameter :: values = 'some message'
   real :: value

   write (*, '(a)')&
   &attr('   <r><W><bo> ERROR: </W>red text on a white background</r>')

   value = 3.4567

!    write (values, '(a)') 'some message'
   write (line, fmt=f) values
   write (*, '(a)') attr(trim(line))

   call alert("info", values, 'value =', value)
   call alert("warn", values, 'value =', value)
   call alert("error", values, 'value =', value)

   z = calculate(5.0, x, y)
   print *, 'x = ', x
   print *, 'y = ', y
   print *, 'z = ', z
contains
   function calculate(a, b, c)
      !! This function change the value of b and c, cause side effect, can't be pure function.
      implicit none
      real, intent(in) :: a
      real, intent(out) :: b, c
      real :: calculate

      b = a*2
      c = a*3
      calculate = a*4
   end function calculate
end program mattr_test
