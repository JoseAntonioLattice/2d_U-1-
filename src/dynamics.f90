module dynamics

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc
  
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)

contains

  subroutine evolution()


  end subroutine evolution
 
  subroutine set_memory(u,L,beta,betai,betaf,nbeta)

    complex(dp), allocatable, dimension(:,:,:) :: u
    integer(i4), intent(in) :: L
    real(dp), allocatable, dimension(:) :: beta
    real(dp), intent(in) :: betai, betaf
    integer(i4), intent(in) :: nbeta

    call set_pbc(L)
    allocate(u(2,L,L))
    allocate(beta(nbeta))

    do i_beta = 1, nbeta 
       beta(i_beta) = betai + (betaf - betai)/(nbeta-1)*(i-1)
    end do
    
  end subroutine set_memory
  
  subroutine initialization(u,beta,N_thermalization,N_measurements, N_skip)

    complex(dp), dimension(:,:,:), intent(out) :: u
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), intent(in) :: beta

    integer(i4) :: i_skip, i_sweeps
    integer(i4) :: L

    L = size(u(1,:,1))
    
    call hot_start(u,L)
    call thermalization(u, N_thermalization, beta)

    do i_sweeps = 1, N_measurements
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       print*, action(u)
    end do
    
  end subroutine initialization
  
  subroutine thermalization(u,N_thermalization,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta

    integer(i4) :: i_sweeps!,L

    !L = size(u(1,:,1))

    !call hot_start(u,L)
    !call cold_start(u)

    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
    end do
    
  end subroutine thermalization

  subroutine hot_start(u,L)

    integer(i4), intent(in) :: L
    complex(dp), intent(out), dimension(2,L,L) :: u
    real(dp), dimension(2,L,L) :: phi

    call random_number(phi)

    phi = 2*pi*phi
    u = exp(i*phi)
    
  end subroutine hot_start

  subroutine cold_start(u)

    complex(dp), intent(out), dimension(:,:,:) :: u

    u = 1.0_dp
    
  end subroutine cold_start
  
  subroutine sweeps(u,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    real(dp), intent(in) :: beta

    integer(i4) :: x,y,mu
    integer(i4) :: L

    L = size(u(1,:,1))
    
    do x = 1, L
       do y = 1, L
          do mu = 1, 2
             call metropolis(U,[x,y],mu,beta)
          end do
       end do
    end do
    
  end subroutine sweeps

  subroutine metropolis(u,x,mu,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4), intent(in) :: x(2), mu
    real(dp), intent(in) :: beta
    
    complex(dp) :: u_new
    real(dp) :: phi
    real(dp) :: p, deltaS
    real(dp) :: r
    
    call random_number(phi)
    phi = 2*pi*phi
    u_new = exp(i*phi)

    deltaS = DS(u(mu,x(1),x(2)),u_new,beta,staples(u,x,mu))

    call random_number(r)
    p = min(1.0_dp,exp(-DeltaS))
    if ( r <= p )then
       u(mu,x(1),x(2)) = u_new
    end if
    
  end subroutine metropolis

  function staples(u,x,mu)

    complex(dp) :: staples
    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4), intent(in) :: x(2), mu
    integer(i4), dimension(2) :: x2, x3, x4, x5, x6

    integer(i4) :: nu
    
    if ( mu == 1 ) then
       nu = 2
    else
       nu = 1
    end if

    x2 = ipf(x,nu)
    x3 = ipf(x,mu)
    x4 = imf(x,nu)
    x5 = x4
    x6 = imf(x3,nu)
    
    staples = u(nu,x(1),x(2)) * u(mu,x2(1), x2(2)) * conjg( u(nu,x3(1), x3(2)) ) + &
         conjg( u(nu,x4(1),x4(2)) ) * u(mu,x5(1), x5(2)) * u(nu,x6(1), x6(2))
    
  end function staples
  
  function DS(uold, unew, beta,stp)

    real(dp) :: DS,beta
    complex(dp) :: uold, unew, stp

    DS = -beta * real( (unew - uold) * conjg(stp) )
    
  end function DS
  
  function plaquette(u,x)

    complex(dp) :: plaquette
    
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2)
    integer(i4), dimension(2) :: x2, x3


    x2 = ipf(x,1)
    x3 = ipf(x,2)
    
    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * &
          conjg(U(1,x3(1),x3(2))) * conjg(U(2,x(1),x(2)))
    
  end function plaquette

  function action(u)

    real(dp) :: action
    
    complex(dp), dimension(:,:,:), intent(in) :: u

    integer(i4) :: x,y
    integer(i4) :: L

    L = size(U(1,:,1))

    action = 0.0_dp
    
    do x = 1, L
       do y = 1, L
          action = action + real(plaquette(u,[x,y]))
       end do
    end do

    action = action / L**2
    
  end function action
  
end module dynamics
