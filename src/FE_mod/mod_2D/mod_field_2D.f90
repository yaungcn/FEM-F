module mod_field_2D
   use quadrature_module, only: wp => quadrature_wp
   implicit none

   private
   public :: field

   type :: field
      !! Field parameter input.
      real(wp) :: left, right, bottom, top
         !! The problem domain is [left,right]*[bottom,top].
      real(wp) :: h_partition, v_partition
         !! The step size of the vertical and horizontal partition.
      integer :: Nh, Nv
         !! The number of horizontal and vertical partition.
      real(wp) :: tolerance = 10000*epsilon(1.0_wp)
      integer :: Gauss_point_number
         !! The number of Gauss Quadrature points.
      integer :: basis_type = 201
         !! 201: 2D linear;
         !! 202: 2D quadratic;
      integer :: mesh_type = 503
         !! 503: 2D triangular mesh;
         !! 504: 2D quadrilateral mesh;

   contains
      procedure, pass(self) :: init => init_field, check => check_field

   end type field

contains
   subroutine init_field(self, left, right, bottom, top, Nh, Nv, Gauss_point_number, basis_type, mesh_type)
      !! TODO: to initialize the field information object 'field'.
      class(field), intent(inout) :: self
      real(wp), intent(in) :: left, right, bottom, top
      integer, intent(in) :: Nh, Nv
      integer, intent(in) :: Gauss_point_number
      integer, intent(in), optional :: basis_type
      integer, intent(in), optional :: mesh_type

      self%left = left
      self%right = right
      self%bottom = bottom
      self%top = top
      self%Nh = Nh
      self%Nv = Nv
      self%h_partition = (right - left)/real(Nh, wp)
      self%v_partition = (top - bottom)/real(Nv, wp)
      self%Gauss_point_number = Gauss_point_number
      select case (basis_type)
      case (201)
         self%basis_type = basis_type
      case (202)
         self%basis_type = basis_type
      case default
         self%basis_type = 201
      end select

      select case (mesh_type)
      case (503)
         self%mesh_type = mesh_type
      case (504)
         self%mesh_type = mesh_type
      case default
         self%basis_type = 503
      end select
   end subroutine init_field

   subroutine check_field(self)
      !! TODO: to check the field information object 'field'.
      class(field), intent(inout) :: self

      if (self%left >= self%right) then
         error stop 'Error: left >= right'

      else if (self%bottom >= self%top) then
         error stop 'Error: bottom >= top'

      else if (self%h_partition <= 0) then
         error stop 'Error: h_partition <= 0'

      else if (self%v_partition <= 0) then
         error stop 'Error: v_partition <= 0'

      else if (self%Nh <= 0) then
         error stop 'Error: Nh <= 0'

      else if (self%Nv <= 0) then
         error stop 'Error: Nv <= 0'

      else if (self%Gauss_point_number <= 0) then
         error stop 'Error: Gauss_point_number <= 0'

      else if (self%basis_type /= 201 .and. self%basis_type /= 202) then
         error stop 'Error: basis_type /= 201(for 2D linear, default) and basis_type /= 202(for 2D quadratic)'

      else if (self%mesh_type /= 503 .and. self%mesh_type /= 504) then
         error stop 'Error: mesh_type /= 503(for triangluar, default) and mesh_type /= 504(for quadrilateral)'

      end if

      write (*, '(A)') 'Field information is checked.'

   end subroutine check_field

end module mod_field_2D
