module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  integer(i4) :: L
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta

  integer(i4) :: inunit,outunit
  character(99) :: inputfilename, outputfilename
    
  
  namelist /parametersfile/ L,N_thermalization,N_measurements,N_skip, &
       beta_i, beta_f, n_beta
contains

  subroutine read_input()

    
    read(*,*) inputfilename
    print*, 'Enter input parameters file: ', trim(inputfilename)
    read(*,*) outputfilename
    print*, 'Enter output parameters file: ', trim(outputfilename)
    open(newunit = inunit,file = trim(inputfilename), status = 'old')
    read(inunit, nml = parametersfile)
    write(*,nml = parametersfile)
  end subroutine read_input

  
end module parameters
