module arrays

  use iso_fortran_env, only : dp => real64

  implicit none

  complex(dp), allocatable, dimension(:,:,:) :: u
  real(dp), allocatable, dimension(:) :: beta 
  
end module arrays
