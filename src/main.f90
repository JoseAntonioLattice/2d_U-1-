program U1_2d

  use iso_fortran_env
  use statistics
  use pbc
  use parameters
  use arrays
  use dynamics
  use number2string_mod
  use create_files
  implicit none

  integer :: i_b
  character(:), allocatable :: datafile
  real(dp) :: r
  integer(int64) :: rate, start_time, end_time
  
  call read_input
  call set_memory(u,L,beta,beta_i,beta_f,n_beta,plq_action,n_measurements)
  print*, beta
  !print*,sqrt(2/beta)
  !call cold_start(u)
  !call random_number(r)
  !call sleep(floor(r*10+1))
  call create_measurements_file_2(L,trim(outputfilename),outunit)
  call system_clock(count_rate = rate)
  call system_clock(count = start_time)
  call hot_start(u,L)
  !open( newunit = outunit, file = outputfilename, status = 'old')
  do i_b = 1, n_beta
     call initialization(u,plq_action,beta(i_b),N_thermalization,N_measurements, N_skip)
     call max_jackknife_error_2(plq_action,avr_action,err_action,bins)
     print*,beta(i_b), avr_action, err_action
     write(outunit,*) beta(i_b), avr_action, err_action
     flush(outunit)
  end do
  call system_clock(count = end_time)
  write(outunit,*) ' '
  write(outunit,*) ' '
  write(outunit,*) ' '
  write(outunit,*) real(end_time-start_time)/real(rate)
  close(outunit)
  

end program U1_2d
