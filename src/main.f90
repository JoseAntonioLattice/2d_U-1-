program U1_2d

  use pbc
  use parameters
  use arrays
  use dynamics
  
  implicit none

  integer :: i_b

  call read_input
  call set_memory(u,L,beta,beta_i,beta_f,n_beta)

  do i_b = 1, n_beta
     call initialization(u,beta(i_b),N_thermalization,N_measurements, N_skip)
  end do

end program U1_2d
