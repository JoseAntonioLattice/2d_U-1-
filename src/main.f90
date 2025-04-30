program U1_2d

  use statistics
  use pbc
  use parameters
  use arrays
  use dynamics
  use number2string_mod
  use create_files
  implicit none

  integer :: i_b, outfile
  character(:), allocatable :: datafile
  real(dp) :: r
  
  call read_input
  call set_memory(u,L,beta,beta_i,beta_f,n_beta,plq_action,n_measurements)
  print*, beta
  !print*,sqrt(2/beta)
  !call cold_start(u)
  call random_number(r)
  !call sleep(floor(r*10+1))
  call create_measurements_file(L,datafile)
  call hot_start(u,L)
  open( newunit = outfile, file = datafile, status = 'old')
  do i_b = 1, n_beta
     call initialization(u,plq_action,beta(i_b),N_thermalization,N_measurements, N_skip)
     call max_jackknife_error_2(plq_action,avr_action,err_action,bins)
     print*,beta(i_b), avr_action, err_action
     write(outfile,*) beta(i_b), avr_action, err_action
     flush(outfile)
  end do
  close(outfile)

end program U1_2d
