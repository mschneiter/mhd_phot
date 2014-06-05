!=======================================================================
!   This module contains global variables
!=======================================================================
module globals
  implicit none
  !
  !   flow variables
  real, dimension(:,:,:,:), allocatable :: u, up, primit
  !   fluxes
  real, dimension(:,:,:,:), allocatable :: f, g, h
  !
  !---------------------------------------------------------------------
  !   spacing
  real :: dx, dy, dz
  !
  !   position of block and neighbors
  integer, dimension(0:2) :: coords
  integer :: left, right, top, bottom, out, in
  !
  !   MPI rank & Cartesian communicator
  integer :: rank, comm2d
  !
end module globals
!=======================================================================
