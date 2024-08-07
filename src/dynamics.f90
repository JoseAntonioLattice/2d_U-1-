module dynamics

  use iso_fortran_env, only : dp => real64, i4 => int32
  
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
  
  subroutine thermalization(u,N_thermalization,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta

    integer(i4) :: i_sweeps

    call hot_start(u)
    !call cold_start(u)
    
  end subroutine thermalization


  subroutine hot_start(u)

    complex(dp), intent(in), dimension(:,:,:) :: u
    integer(i4), parameter :: d = size(u(:,1,1))
    integer(i4), parameter :: L = size(u(1,:,1))
    real(dp), dimension(d,L,L) :: phi

    call random_number(phi)

    phi = 2*pi*phi
    u = exp(i*phi)
    
  end subroutine hot_start

  subroutine cold_start(u)

    complex(dp), intent(in), dimension(:,:,:) :: u

    u = 1.0_dp
    
  end subroutine cold_start
  
  subroutine sweeps(u,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    real(dp), intent(in) :: beta

    integer(i4) :: x,y,mu

    integer(i4), parameter :: d = size(u(:,1,1))
    integer(i4), parameter :: L = size(u(1,:,1))
    
    do x = 1, L
       do y = 1, L
          do mu = 1, d
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

    deltaS = DS(u(mu,x(1),x(2)),u_new,staples(u,x,mu))

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
    integer(i4), dimension(2), intent(in) :: x2, x3, x4, x5, x6
    
    if ( mu == 1 ) then
       nu = 2
    else
       nu = 1
    end if
    
    staples = u(nu,x(1),x(2)) * u(mu,x2(1), x2(2)) * conjg( u(nu,x3(1), x3(2)) ) + &
         conjg( u(nu,x4(1),x4(2)) ) * u(mu,x5(1), x5(2)) * u(nu,x6(1), x6(2))
    
  end function staples
  
  function DS(uold, unew, stp)

    real(dp) :: DS
    complex(dp) :: uold, unew, stp

    DS = -beta * real( (unew - uold) * conjg(stp) )
    
  end function DS
  
  function plaquette(u,x)

    complex(dp) :: plaquette
    
    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4), intent(in) :: x(2), mu

    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * conjg(U(1,x3(1),x3(2))) * conjg(U(2,x4(1),x4(2)))
    
  end function plaquette
  
end module dynamics
