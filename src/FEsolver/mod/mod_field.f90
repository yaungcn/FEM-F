module mod_field
   use quadrature_module, only: wp => quadrature_wp
   implicit none

   private
   public :: field
   
   type :: field
      !! Field parameter input.
      real(wp) :: left, right
         !! The problem domain is [left,right]*[bottom,top].
      real(wp) :: h_partition
         !! The step size of the partition.
      integer :: Gauss_point_number
         !! The number of Gauss Quadrature points.
      integer :: basis_type = 101
         !! 101: 1D linear;
         !! 102: ?? ???

   contains
      procedure, pass(self) :: init => init_field

   end type field

contains
   subroutine init_field(self, left, right, h_partition, Gauss_point_number, basis_type)
      !! TODO: to initialize the field information object 'field'.
      class(field), intent(inout) :: self
      real(wp), intent(in) :: left, right, h_partition
      integer, intent(in) :: Gauss_point_number
      integer, intent(in), optional :: basis_type

      self%left = left
      self%right = right
      self%h_partition = h_partition
      self%Gauss_point_number = Gauss_point_number
      select case (basis_type)
      case (101)
         self%basis_type = basis_type
      case (102)
         self%basis_type = basis_type
      case default
         self%basis_type = 101
      end select
   end subroutine init_field

end module mod_field
