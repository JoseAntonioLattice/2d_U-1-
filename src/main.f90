program U1_2d

  use pbc
  use parameters
  use arrays
  use dynamics
  
  implicit none

  call read_input
  call set_memory(u,L)
  call initialization(u,beta_f,N_thermalization,N_measurements, N_skip)

end program U1_2d
