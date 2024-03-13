module mod_field_2D
   use mod_kinds
   use mod_message, only: info_print, warn_print, error_print
   use M_attr, only: attr
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
         !! The tolerance for the iterative solver.
      integer :: Gauss_point_number = 4
         !! The number of Gauss Quadrature points.
      integer :: trial_basis_type = 201
         !! 201: 2D linear;
         !! 202: 2D quadratic;
      integer :: test_basis_type = 201
         !! 201: 2D linear;
         !! 202: 2D quadratic;
      integer :: mesh_type = 503
         !! 503: 2D triangular mesh;
         !! 504: 2D quadrilateral mesh;
      integer :: number_of_elements = 0
         !! The number of elements in the mesh.
      integer :: number_of_nodes_mesh = 0
         !! The number of nodes in the mesh.
      integer :: number_of_nodes_fe = 0
         !! The number of nodes in the FE mesh.
      integer :: number_of_local_basis_trial = 0
         !! The number of local basis functions for trial function.
      integer :: number_of_local_basis_test = 0
         !! The number of local basis functions for test function.

   contains
      procedure, pass(self) :: init => init_field, check => check_field, print => print_field

   end type field

contains
   subroutine init_field(self, left, right, bottom, top, &
                         Nh_partition, Nv_partition, &
                         Nh_basis, Nv_basis, &
                         Gauss_point_number, tolerance, &
                         trial_basis_type, test_basis_type, mesh_type)
      !! TODO: to initialize the field information object 'field'.
      class(field), intent(inout) :: self
      real(wp), intent(in) :: left, right, bottom, top
      integer, intent(in) :: Nh_partition, Nv_partition
      integer, intent(in) :: Nh_basis, Nv_basis
      integer, intent(in), optional :: Gauss_point_number
      real(wp), intent(in), optional :: tolerance
      integer, intent(in), optional :: trial_basis_type, test_basis_type
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

      self%number_of_nodes_mesh = (Nh_partition + 1)*(Nv_partition + 1)
      self%number_of_nodes_fe = (Nh_basis + 1)*(Nv_basis + 1)

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

      if (present(trial_basis_type)) then
         self%trial_basis_type = trial_basis_type
      else
         self%trial_basis_type = 201
         call warn_print('trial_basis_type not set, use default value 201(for 2D linear) or 202(for 2D quadratic)')
         ! write (*, '(A)') 'Warning: basis_type /= 201(for 2D linear) or 202(for 2D quadratic), use default value 201'
      end if

      if (present(test_basis_type)) then
         self%test_basis_type = test_basis_type
      else
         self%test_basis_type = 201
         call warn_print('test_basis_type not set, use default value 201(for 2D linear) or 202(for 2D quadratic)')
         ! write (*, '(A)') 'Warning: basis_type /= 201(for 2D linear) or 202(for 2D quadratic), use default value 201'
      end if

      if (present(mesh_type)) then
         self%mesh_type = mesh_type
      else
         self%mesh_type = 503
         call warn_print('mesh_type not set, use default value 503(for triangluar) or 504(for quadrilateral)')
         ! write (*, '(A)') 'Warning: mesh_type /= 503(for triangluar) or 504(for quadrilateral), use default value 503'
      end if

      select case (mesh_type)
      case (503)
         self%number_of_elements = 2*Nh_partition*Nv_partition
         select case (trial_basis_type)
         case (201)
            self%number_of_local_basis_trial = 3
         case (202)
            self%number_of_local_basis_trial = 6
         case default
            call error_print('Undefined trial_basis_type')
            stop
         end select
         select case (test_basis_type)
         case (201)
            self%number_of_local_basis_test = 3
         case (202)
            self%number_of_local_basis_test = 6
         case default
            call error_print('Undefined test_basis_type')
            stop
         end select
      case (504)
         self%number_of_elements = Nh_partition*Nv_partition
         select case (trial_basis_type)
         case (201)
            self%number_of_local_basis_trial = 4
         case (202)
            self%number_of_local_basis_trial = 9
         case default
            call error_print('Undefined trial_basis_type')
            stop
         end select
         select case (test_basis_type)
         case (201)
            self%number_of_local_basis_test = 4
         case (202)
            self%number_of_local_basis_test = 9
         case default
            call error_print('Undefined test_basis_type')
            stop
         end select
      case default
         call error_print('Undefined mesh_type')
         stop
      end select

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

      else if (self%trial_basis_type /= 201 .and. self%trial_basis_type /= 202) then
         error stop 'Error: trial_basis_type /= 201(for 2D linear, default) and basis_type /= 202(for 2D quadratic)'

      else if (self%test_basis_type /= 201 .and. self%test_basis_type /= 202) then
         error stop 'Error: test_basis_type /= 201(for 2D linear, default) and basis_type /= 202(for 2D quadratic)'

      else if (self%mesh_type /= 503 .and. self%mesh_type /= 504) then
         error stop 'Error: mesh_type /= 503(for triangluar, default) and mesh_type /= 504(for quadrilateral)'

      else if (self%number_of_elements <= 0) then
         error stop 'Error: number_of_elements <= 0'
      else if (self%number_of_nodes_mesh <= 0) then
         error stop 'Error: number_of_nodes_mesh <= 0'
      else if (self%number_of_nodes_fe <= 0) then
         error stop 'Error: number_of_nodes_fe <= 0'
      else if (self%number_of_local_basis_trial <= 0) then
         error stop 'Error: number_of_local_basis_trial <= 0'
      else if (self%number_of_local_basis_test <= 0) then
         error stop 'Error: number_of_local_basis_test <= 0'

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
      write (*, '(A)') attr(' |<bo><in> Field information: </in></bo>')
      write (*, float_format) ' |   left, right = ', self%left, self%right
      write (*, float_format) ' |   bottom, top = ', self%bottom, self%top
      write (*, '(A, 2I6)') ' |   Nh_partition, Nv_partition = ', self%Nh_partition, self%Nv_partition
      write (*, float_format) ' |   ->h_partition, v_partition = ', self%h_partition, self%v_partition
      write (*, '(A, 2I6)') ' |   Nh_basis, Nv_basis = ', self%Nh_basis, self%Nv_basis
      write (*, float_format) ' |   ->h_basis, v_basis = ', self%h_basis, self%v_basis
      write (*, '(A, I6)') ' |   Gauss_point_number = ', self%Gauss_point_number
      write (*, '(A, 2F12.6)') ' |   tolerance  = ', self%tolerance
      write (*, '(A, I6)') ' |   trial_basis_type = ', self%trial_basis_type
      write (*, '(A, I6)') ' |   test_basis_type  = ', self%test_basis_type
      write (*, '(A, I6)') ' |   mesh_type  = ', self%mesh_type
      write (*, '(A, I6)') ' |   number_of_elements   = ', self%number_of_elements
      write (*, '(A, I6)') ' |   number_of_nodes_mesh = ', self%number_of_nodes_mesh
      write (*, '(A, I6)') ' |   number_of_nodes_fe   = ', self%number_of_nodes_fe
      write (*, '(A, I6)') ' |   number_of_local_basis_trial = ', self%number_of_local_basis_trial
      write (*, '(A, I6)') ' |   number_of_local_basis_test  = ', self%number_of_local_basis_test
      write (*, '(A)') ' -----------------------------------------------------------------'

   end subroutine print_field
end module mod_field_2D
