program U1_2d

  use statistics
  use pbc
  use parameters
  use arrays
  use dynamics
  
  implicit none

  integer :: i_b

  call read_input
  call set_memory(u,L,beta,beta_i,beta_f,n_beta,plq_action,n_measurements)
  print*, beta
  !print*,sqrt(2/beta)
  call hot_start(u,L)
  open( unit = 10, file = 'data/data.dat', status = 'unknown')
  do i_b = 1, n_beta
     call initialization(u,plq_action,beta(i_b),N_thermalization,N_measurements, N_skip)
     call max_jackknife_error_2(plq_action,avr_action,err_action,bins)
     print*,beta(i_b), avr_action, err_action
     write(10,*) beta(i_b), avr_action, err_action
     flush(10)
  end do

end program U1_2d
