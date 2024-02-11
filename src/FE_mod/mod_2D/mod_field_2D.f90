module mod_field_2D
   use quadrature_module, only: wp => quadrature_wp
   use mod_message, only: info_print, warn_print, error_print
   implicit none

   private
   public :: field

   type :: field
      !! Field parameter input.
      real(wp) :: left, right, bottom, top
         !! The problem domain is [left,right]*[bottom,top].
      real(wp) :: h_partition, v_partition
         !! The step size of the vertical and horizontal partition.
         !! h_partition = (right - left)/Nh
         !! v_partition = (top - bottom)/Nv.
      real(wp) :: h_basis, v_basis
         !! The step size of the vertical and horizontal FE basis func.
         !! h_basis = (right - left)/Nh_basis
         !! v_basis = (top - bottom)/Nv_basis.
      integer :: Nh_partition, Nv_partition
         !! The number of horizontal and vertical partition.
      integer :: Nh_basis, Nv_basis
         !! The number of the FE basis functions in the horizontal and vertical direction.
      real(wp) :: tolerance = 10000*epsilon(1.0_wp)
         !! The tolerance for the iterative solver.                        Not used yet.!!
      integer :: Gauss_point_number = 4
         !! The number of Gauss Quadrature points.
      integer :: basis_type = 201
         !! 201: 2D linear;
         !! 202: 2D quadratic;
      integer :: mesh_type = 503
         !! 503: 2D triangular mesh;
         !! 504: 2D quadrilateral mesh;

   contains
      procedure, pass(self) :: init => init_field, check => check_field, print => print_field

   end type field

contains
   subroutine init_field(self, left, right, bottom, top, &
                         Nh_partition, Nv_partition, &
                         Nh_basis, Nv_basis, &
                         Gauss_point_number, tolerance, basis_type, mesh_type)
      !! TODO: to initialize the field information object 'field'.
      class(field), intent(inout) :: self
      real(wp), intent(in) :: left, right, bottom, top
      integer, intent(in) :: Nh_partition, Nv_partition
      integer, intent(in) :: Nh_basis, Nv_basis
      integer, intent(in), optional :: Gauss_point_number
      real(wp), intent(in), optional :: tolerance
      integer, intent(in), optional :: basis_type
      integer, intent(in), optional :: mesh_type

      self%left = left
      self%right = right
      self%bottom = bottom
      self%top = top
      self%Nh_partition = Nh_partition
      self%Nv_partition = Nv_partition
      self%Nh_basis = Nh_basis
      self%Nv_basis = Nv_basis
      self%h_partition = (right - left)/real(Nh_partition, wp)
      self%v_partition = (top - bottom)/real(Nv_partition, wp)
      self%h_basis = (right - left)/real(Nh_basis, wp)
      self%v_basis = (top - bottom)/real(Nv_basis, wp)

      if (present(Gauss_point_number)) then
         self%Gauss_point_number = Gauss_point_number
      else
         self%Gauss_point_number = 4
         call warn_print('Gauss_point_number not set, use default value 4')
         ! write (*, '(3x,a)') attr('<B><bo><y> WARNING </B></y>  Gauss_point_number not set, use default value 4')
         ! write (*, '(A)') 'Warning: Gauss_point_number not set, use default value 4'
      end if

      if (present(tolerance)) then
         self%tolerance = tolerance
      else
         self%tolerance = 1.0e-2_wp
         call warn_print('tolerance not set, use default value 1.0e-2')
         ! write (*, '(3x,a)') attr('<B><bo><y> WARNING </B></y>  tolerance not set, use default value 1.0e-2')
         ! write (*, '(A)') 'Warning: tolerance not set, use default value 1.0e-2'
      end if

      if (present(basis_type)) then
         self%basis_type = basis_type
      else
         self%basis_type = 201
         call warn_print('basis_type not set, use default value 201(for 2D linear) or 202(for 2D quadratic)')
         ! write (*, '(A)') 'Warning: basis_type /= 201(for 2D linear) or 202(for 2D quadratic), use default value 201'
      end if

      if (present(mesh_type)) then
         self%mesh_type = mesh_type
      else
         self%mesh_type = 503
         call warn_print('mesh_type not set, use default value 503(for triangluar) or 504(for quadrilateral)')
         ! write (*, '(A)') 'Warning: mesh_type /= 503(for triangluar) or 504(for quadrilateral), use default value 503'
      end if

      call info_print('Field information initialized.')
      call check_field(self)

   end subroutine init_field

   subroutine check_field(self)
      !> TODO: to check the field information object 'field'.
      class(field), intent(inout) :: self

      if (self%left >= self%right) then
         error stop 'Error: left >= right'

      else if (self%bottom >= self%top) then
         error stop 'Error: bottom >= top'

      else if (self%h_partition <= 0) then
         error stop 'Error: h_partition <= 0'

      else if (self%v_partition <= 0) then
         error stop 'Error: v_partition <= 0'

      else if (self%h_basis <= 0) then
         error stop 'Error: h_basis <= 0'

      else if (self%v_basis <= 0) then
         error stop 'Error: v_basis <= 0'

      else if (self%Nh_partition <= 0) then
         error stop 'Error: Nh <= 0'

      else if (self%Nv_partition <= 0) then
         error stop 'Error: Nv <= 0'

      else if (self%Nh_basis <= 0) then
         error stop 'Error: Nh_basis <= 0'

      else if (self%Nv_basis <= 0) then
         error stop 'Error: Nv_basis <= 0'

      else if (self%tolerance <= 0) then
         error stop 'Error: tolerance <= 0'

      else if (self%Gauss_point_number <= 0) then
         error stop 'Error: Gauss_point_number <= 0'

      else if (self%basis_type /= 201 .and. self%basis_type /= 202) then
         error stop 'Error: basis_type /= 201(for 2D linear, default) and basis_type /= 202(for 2D quadratic)'

      else if (self%mesh_type /= 503 .and. self%mesh_type /= 504) then
         error stop 'Error: mesh_type /= 503(for triangluar, default) and mesh_type /= 504(for quadrilateral)'

      end if

      call info_print('Field information checked.')

   end subroutine check_field

   subroutine print_field(self, float_format_string)
      !> TODO: to print the field information object 'field'.
      class(field), intent(in) :: self
      character(len=*), intent(in), optional :: float_format_string
      !! DEFAULT: float_format = '(A, 2F12.6)'
      character(len=100) :: float_format

      if (present(float_format_string)) then
         float_format = float_format_string
      else
         float_format = '(A, 2F12.6)'
      end if

      write (*, '(A)') ' -----------------------------------------------------------------'
      write (*, '(A)') ' | Field information:'
      write (*, float_format) ' |   left, right = ', self%left, self%right
      write (*, float_format) ' |   bottom, top = ', self%bottom, self%top
      write (*, '(A, 2I6)') ' |   Nh_partition, Nv_partition = ', self%Nh_partition, self%Nv_partition
      write (*, float_format) ' |   ->h_partition, v_partition = ', self%h_partition, self%v_partition
      write (*, '(A, 2I6)') ' |   Nh_basis, Nv_basis = ', self%Nh_basis, self%Nv_basis
      write (*, float_format) ' |   ->h_basis, v_basis = ', self%h_basis, self%v_basis
      write (*, '(A, I6)') ' |   Gauss_point_number = ', self%Gauss_point_number
      write (*, '(A, 2F12.6)') ' |   tolerance  = ', self%tolerance
      write (*, '(A, I6)') ' |   basis_type = ', self%basis_type
      write (*, '(A, I6)') ' |   mesh_type  = ', self%mesh_type
      write (*, '(A)') ' -----------------------------------------------------------------'

   end subroutine print_field
end module mod_field_2D
