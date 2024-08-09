module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  integer(i4) :: L
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta
  
  namelist /parametersfile/ L,N_thermalization,N_measurements,N_skip, &
       beta_i, beta_f, n_beta
contains

  subroutine read_input()

    integer(i4) :: inunit
    character(99) :: inputfilename
    
    read(*,*) inputfilename
    print*, 'User typed: ', trim(inputfilename)
    open(newunit = inunit,file = trim(inputfilename), status = 'old')
    read(inunit, nml = parametersfile)
    write(*,nml = parametersfile)
  end subroutine read_input

  
end module parameters
