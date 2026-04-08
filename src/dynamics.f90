module dynamics

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc
  use parameters, only : L
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)

contains
 
  subroutine set_memory(u,beta,betai,betaf,nbeta, plqaction, corr_poly, n_measurements)

    complex(dp), allocatable, dimension(:,:,:) :: u
    real(dp), allocatable, dimension(:) :: beta
    real(dp), intent(in) :: betai, betaf
    real(dp), allocatable, dimension(:) :: plqaction
    integer(i4), intent(in) :: nbeta, n_measurements
    real(dp), allocatable, dimension(:) :: beta_copy
    integer(i4) :: i_beta
    complex(dp), allocatable :: corr_poly(:,:)
    
    call set_pbc(L)
    allocate(u(2,L(1),L(2)))
    allocate(beta(nbeta), beta_copy(nbeta))

    do i_beta = 1, nbeta 
       beta(i_beta) = betai + (betaf - betai)/(nbeta-1)*(i_beta-1)
    end do
    !beta = 1/beta**2
    !beta_copy = beta

   !do i_beta = 1, nbeta
   !    beta(i_beta) = beta_copy(nbeta+1-i_beta)
   !end do
    
    allocate(plqaction(n_measurements))
    allocate(corr_poly(0:L(1)/2-1,n_measurements))
    
  end subroutine set_memory
  
  subroutine initialization(u,plqaction,corr_poly,beta,N_thermalization,N_measurements, N_skip)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    real(dp), dimension(:), intent(out) :: plqaction
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), intent(in) :: beta
    complex(dp) :: corr_poly(0:L(1)/2-1,N_measurements)
    integer(i4) :: i_skip, i_sweeps

    call thermalization(u, N_thermalization, beta)

    do i_sweeps = 1, N_measurements
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       plqaction(i_sweeps) = action(u)
       !corr_poly(:,i_sweeps) = correlation_polyakov(U)
       corr_poly(0:,i_sweeps) = correlation_polyakov(U)
    end do
    
  end subroutine initialization
  
  subroutine thermalization(u,N_thermalization,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta

    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
    end do
    
  end subroutine thermalization

  subroutine hot_start(u)

    complex(dp), intent(out), dimension(2,L(1),L(2)) :: u
    real(dp), dimension(2,L(1),L(2)) :: phi

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
    
    !call hmc(U,beta,50,1.0_dp,L)
    do x = 1, L(1)
       do y = 1, L(2)
          do mu = 1, 2
             call heatbath(U,[x,y],mu,beta)
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

  subroutine heatbath(u,x,mu,beta)
     complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4), intent(in) :: x(2), mu
    real(dp), intent(in) :: beta
    
    complex(dp) :: sigma
    real(dp) :: phi, s, gamma
    real(dp) :: rho_max
    real(dp) :: r
    logical :: condition
    
    sigma = staples(u,x,mu)
    s = abs(sigma)
    gamma = atan2(sigma%im,sigma%re)
    rho_max = exp(beta*s)

    condition = .false.
    do while(.not.condition)
       call random_number(phi)
       phi = (2*phi - 1.0_dp)*pi
       call random_number(r)
       r = rho_max*r
       if( r <= rho_max**(cos(phi-gamma)) ) condition = .true.
    end do

    U(mu,x(1),x(2)) = exp(i*phi)
    
  end subroutine heatbath
  
  
  subroutine hmc(U, beta, N, Time)
    complex(dp), intent(inout) :: U(2,L(1),L(2))
    integer(i4), intent(in) :: N
    real(dp), intent(in) :: beta, Time
    complex(dp), dimension(2,L(1),L(2)) ::  Up
    real(dp), dimension(2,L(1),L(2)) :: Forces, p, pnew
    integer(i4) :: k, x, y, mu
    real(dp) :: DeltaH,r, deltaT

    deltat = Time/N
    call generate_pi(p)

    !! k = 0
    !U_0
    up = u
    !P_0
    pnew = p

    !Compute F[U_0]
    call compute_forces(Forces,beta,u)

    ! Compute P_{1/2} = P_0 + 0.5*dt*F[U_0]
    pnew = pnew + 0.5*deltaT*Forces
 
    ! k = 1, n -1
    do k = 1, N - 1
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       up = up * exp(i*DeltaT*pnew)
       
       !compute F[U_k]
       call compute_forces(Forces,beta,up)
      
       !P_{k+1/2} = P_{k-1/2} + dt*F[U_k]
       pnew = pnew + deltaT*Forces
    end do 

    ! k = n
    !U_n = exp(i*dt*P_{n-1/2})U_{n-1}
    up = up * exp(i*DeltaT*pnew)

    !compute F[U_n]
    call compute_forces(Forces,beta,up)
    
    !P_n = P_{n-1/2} + 0.5*dt*F[U_n]
    p = p + 0.5*DeltaT*Forces

    ! Metropolis step
    DeltaH = DH(u,up,p,pnew,beta)
    call random_number(r)
    if( r <= exp(-DeltaH) ) u = up
    
  end subroutine hmc

 
  
  function DH(U,Unew,P,Pnew,beta)
    real(dp) :: DH
    complex(dp), dimension(:,:,:), intent(in) :: U, Unew
    real(dp), dimension(:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x, y,mu
    real(dp) :: DeltaS

    DH = 0.0_dp
    DeltaS = 0.0_dp
    do x = 1, L(1)
       do y = 1, L(2)
          DeltaS = DeltaS + real(plaquette(u,[x,y]) - plaquette(unew,[x,y]))
       end do
    end do

    DeltaS = beta*DeltaS

    do x = 1, L(1)
       do y = 1, L(2)
          do mu = 1, 2
             DH = DH + (pnew(mu,x,y))**2 - (p(mu,x,y))**2
          end do
       end do
    end do

    DH = 0.5*DH + DeltaS
    
    
  end function DH

  
  subroutine generate_pi(p)
    real(dp), intent(out), dimension(2,L(1),L(2)) :: p
    real(dp) :: u1, u2
    integer(i4) :: x, y

    do x = 1, L(1)
       do y = 1, L(2)
          call random_number(u1)
          call random_number(u2)
          p(1,x,y) = sqrt(-2*log(u1))*cos(2*pi*u2) 
          p(2,x,y) = sqrt(-2*log(u1))*sin(2*pi*u2)
       end do
    end do
  end subroutine generate_pi

  subroutine compute_forces(forces,beta,u)
    complex(dp), intent(in) :: u(2,L(1),L(2))
    real(dp), intent(out) :: Forces(2,L(1),L(2))
    real(dp), intent(in) :: beta
    complex(dp) :: stp
    integer(i4) :: x, y, mu

    do x = 1, L(1)
       do y = 1, L(2)
          do mu = 1, 2
             stp = u(mu,x,y)*conjg(staples(u,[x,y],mu))
             forces(mu,x,y) = -beta*stp%im
          end do
       end do
    end do
    
  end subroutine compute_forces
  
  function staples(u,x,mu)

    complex(dp) :: staples
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2), mu
    integer(i4), dimension(2) :: x2, x3, x4, x5, x6

    integer(i4) :: nu
    
    if ( mu == 1 ) then
       nu = 2
    elseif( mu == 2)then
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


    action = 0.0_dp
    
    do x = 1, L(1)
       do y = 1, L(2)
          action = action + real(plaquette(u,[x,y]))
       end do
    end do

    action = action / product(L)
    
  end function action

  
  function action2(u,beta)

    real(dp) :: action2
    
    complex(dp), dimension(:,:,:), intent(in) :: u
    real(dp), intent(in) :: beta
    integer(i4) :: x,y
   
    
    action2 = 0.0_dp
    
    do x = 1, L(1)
       do y = 1, L(2)
          action2 = action2 + beta*(1.0_dp - real(plaquette(u,[x,y])))
       end do
    end do

    !action2 = beta*(L**2 - action2)
    
  end function action2

  function correlation_polyakov(U) result(corr_poly)
    complex(dp), intent(in) :: U(2,L(1),L(2))
    integer(i4) :: x, r, t
    complex(dp), dimension(0:L(1)/2-1) :: corr_poly
    complex(dp), dimension(L(1)) :: polyakov_loop

    do x = 1, L(1)
       polyakov_loop(x) = product(U(2,x,:))
    end do

    corr_poly = 0.0_dp
    do t = 0, L(1)/2 - 1
       do x = 1, L(1)
          r = mod(x-1+t,L(1))+1
          corr_poly(t) = corr_poly(t)+polyakov_loop(x)*conjg(polyakov_loop(r))
       end do
    end do
    corr_poly = corr_poly/L(1)
  end function correlation_polyakov


!=======================================================================
!  FUNCIÓN: I0(x) — Bessel modificada orden 0 (serie de Taylor)
!=======================================================================
real(8) function bessel_i0(x)
  implicit none
  real(8), intent(in) :: x
  real(8) :: term, hx
  integer :: k
  hx = 0.5d0*x;  term = 1.0d0;  bessel_i0 = 1.0d0
  do k = 1, 50
    term = term * (hx/dble(k))**2
    bessel_i0 = bessel_i0 + term
    if (abs(term) < 1.0d-15*abs(bessel_i0)) exit
  end do
end function bessel_i0


!=======================================================================
!  FUNCIÓN: I1(x) — Bessel modificada orden 1 (serie de Taylor)
!=======================================================================
real(8) function bessel_i1(x)
  implicit none
  real(8), intent(in) :: x
  real(8) :: term, hx
  integer :: k
  hx = 0.5d0*x;  term = hx;  bessel_i1 = hx
  do k = 1, 50
    term = term * hx**2 / (dble(k)*dble(k+1))
    bessel_i1 = bessel_i1 + term
    if (abs(term) < 1.0d-15*abs(bessel_i1)) exit
  end do
end function bessel_i1


end module dynamics
