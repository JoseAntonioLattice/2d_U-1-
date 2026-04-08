module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  complex(dp), allocatable, dimension(:,:,:) :: u
  real(dp), allocatable, dimension(:) :: beta 
  real(dp), allocatable, dimension(:) :: plq_action
  complex(dp), allocatable :: corr_poly(:,:)
  real(dp) :: avr_action,err_action
  integer(i4) :: bins
  
end module arrays
