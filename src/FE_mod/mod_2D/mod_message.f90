module mod_message
   use M_attr, only: attr, attr_mode, alert
   implicit none
   private

   public :: info_print, warn_print, error_print

contains

   subroutine info_print(message)
      character(len=*) :: message
      call alert("info", message)
   end subroutine info_print

   subroutine warn_print(message)
      character(len=*) :: message
      call alert("warn", message)
   end subroutine warn_print

   subroutine error_print(message)
      character(len=*) :: message
      call alert("error", message)
   end subroutine error_print

end module mod_message
