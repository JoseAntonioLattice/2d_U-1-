module dynamics

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc
  
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)

contains
 
  subroutine set_memory(u,L,beta,betai,betaf,nbeta, plqaction,n_measurements)

    complex(dp), allocatable, dimension(:,:,:) :: u
    integer(i4), intent(in) :: L
    real(dp), allocatable, dimension(:) :: beta
    real(dp), intent(in) :: betai, betaf
    real(dp), allocatable, dimension(:) :: plqaction
    integer(i4), intent(in) :: nbeta, n_measurements
    real(dp), allocatable, dimension(:) :: beta_copy
    integer(i4) :: i_beta
    
    call set_pbc(L)
    allocate(u(2,L,L))
    allocate(beta(nbeta), beta_copy(nbeta))

    do i_beta = 1, nbeta 
       beta(i_beta) = betai + (betaf - betai)/(nbeta-1)*(i_beta-1)
    end do
    !beta = 2/beta**2
    !beta_copy = beta

    !do i_beta = 1, nbeta
    !   beta(i_beta) = beta_copy(nbeta+1-i_beta)
    !end do
    
    allocate(plqaction(n_measurements))
    
  end subroutine set_memory
  
  subroutine initialization(u,plqaction,beta,N_thermalization,N_measurements, N_skip)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    real(dp), dimension(:), intent(out) :: plqaction
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), intent(in) :: beta

    integer(i4) :: i_skip, i_sweeps

    call thermalization(u, N_thermalization, beta)

    do i_sweeps = 1, N_measurements
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       plqaction(i_sweeps) = action(u)
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

    !call hmc(U,beta,50,1.0_dp,L)
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

  subroutine hmc(U, beta, N, Time, L)
    complex(dp), intent(inout) :: U(2,L,L)
    integer(i4), intent(in) :: N, L
    real(dp), intent(in) :: beta, Time
    complex(dp), dimension(2,L,L) ::  Up
    real(dp), dimension(2,L,L) :: Forces, p, pnew
    integer(i4) :: k, x, y, mu
    real(dp) :: DeltaH,r, deltaT

    deltat = Time/N
    call generate_pi(p,L)

    !! k = 0
    !U_0
    up = u
    !P_0
    pnew = p

    !Compute F[U_0]
    call compute_forces(Forces,beta,u,L)

    ! Compute P_{1/2} = P_0 + 0.5*dt*F[U_0]
    pnew = pnew + 0.5*deltaT*Forces
 
    ! k = 1, n -1
    do k = 1, N - 1
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       up = up * exp(i*DeltaT*pnew)
       
       !compute F[U_k]
       call compute_forces(Forces,beta,up,L)
      
       !P_{k+1/2} = P_{k-1/2} + dt*F[U_k]
       pnew = pnew + deltaT*Forces
    end do 

    ! k = n
    !U_n = exp(i*dt*P_{n-1/2})U_{n-1}
    up = up * exp(i*DeltaT*pnew)

    !compute F[U_n]
    call compute_forces(Forces,beta,up,L)
    
    !P_n = P_{n-1/2} + 0.5*dt*F[U_n]
    p = p + 0.5*DeltaT*Forces

    ! Metropolis step
    DeltaH = DH(u,up,p,pnew,beta)
    call random_number(r)
    if( r <= exp(-DeltaH) ) u = up
    
  end subroutine hmc

  subroutine hmc_jaime(U,beta,Ntime,time,L)
    integer(i4), intent(in) :: L, Ntime
    complex(dp), dimension(2,L,L), intent(inout) :: U
    real(dp), intent(in) :: beta, time

    complex(dp), dimension(2,L,L) :: unew
    real(dp), dimension(2,L,L) :: p, pnew, force

    real(dp) :: r,u1,u2, DeltaH, dt
    integer(i4) :: x, y, mu, k


    dt = time/NTime
    call generate_pi(p,L)
     
    unew = u
    pnew = p

    unew = unew*exp(0.5*i*dt*pnew)
    call compute_forces(force,beta,unew,L)
    
    do k = 1, Ntime - 2
       pnew = pnew + dt*force
       unew = unew*exp(i*dt*pnew)
       call compute_forces(force,beta,unew,L)
    end do
    
    pnew = pnew + dt*force
    unew = unew*exp(0.5*i*dt*pnew)
    
    call random_number(r)
    
    DeltaH = DH(u,unew,p,pnew,beta)
    if( r <= exp(-DeltaH)) u = unew
    
  end subroutine hmc_jaime

  subroutine hmc_knechtli(u,beta,nsteps,time,L)
    integer(i4), intent(in) :: L, nsteps
    complex(dp), intent(inout), dimension(2,L,L) :: u
    real(dp), intent(in) :: time, beta
    complex(dp), dimension(2,L,L) :: unew
    real(dp), dimension(2,L,L) :: p, pnew, force 
    real(dp) :: r, DeltaH, dt
    integer(i4) :: x, y, mu, k

    dt = time/nsteps
    
    call generate_pi(p,L)
    
    unew = u
    pnew = p

    call compute_forces(force,beta,u,L)
    do k = 1, nsteps
       !P_{k-1/2} = P_{k-1} +F_{k-1}        
       pnew = pnew + 0.5*dt*Force
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       unew = unew(mu,x,y) * exp(dt*i*pnew)
       !F_k
       call compute_forces(force,beta,unew,L)
       ! P_k = P_{k-1/2} + dt*F_k
       pnew = pnew + 0.5*dt*Force
    end do

    !Metropolis step
    call random_number(r)
    if( r <= exp(-DeltaH)) u = unew
    
  end subroutine hmc_knechtli

  
  function DH(U,Unew,P,Pnew,beta)
    real(dp) :: DH
    complex(dp), dimension(:,:,:), intent(in) :: U, Unew
    real(dp), dimension(:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x, y,mu, L
    real(dp) :: DeltaS
    L = size(U(1,:,1))
    DH = 0.0_dp
    DeltaS = 0.0_dp
    do x = 1, L
       do y = 1, L
          DeltaS = DeltaS + real(plaquette(u,[x,y]) - plaquette(unew,[x,y]))
       end do
    end do

    DeltaS = beta*DeltaS

    do x = 1, L
       do y = 1, L
          do mu = 1, 2
             DH = DH + (pnew(mu,x,y))**2 - (p(mu,x,y))**2
          end do
       end do
    end do

    DH = 0.5*DH + DeltaS
    
    
  end function DH

  
  subroutine generate_pi(p,L)
    integer(i4), intent(in) :: L
    real(dp), intent(out), dimension(2,L,L) :: p
    real(dp) :: u1, u2
    integer(i4) :: x, y

    do x = 1, L
       do y = 1, L
          call random_number(u1)
          call random_number(u2)
          p(1,x,y) = sqrt(-2*log(u1))*cos(2*pi*u2) 
          p(2,x,y) = sqrt(-2*log(u1))*sin(2*pi*u2)
       end do
    end do
  end subroutine generate_pi

  subroutine compute_forces(forces,beta,u,L)
    integer(i4), intent(in) :: L
    complex(dp), intent(in) :: u(2,L,L)
    real(dp), intent(out) :: Forces(2,L,L)
    real(dp), intent(in) :: beta
    complex(dp) :: stp
    integer(i4) :: x, y, mu

    do x = 1, L
       do y = 1, L
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

  
  function action2(u,beta)

    real(dp) :: action2
    
    complex(dp), dimension(:,:,:), intent(in) :: u
    real(dp), intent(in) :: beta
    integer(i4) :: x,y
    integer(i4) :: L

    L = size(U(1,:,1))

    action2 = 0.0_dp
    
    do x = 1, L
       do y = 1, L
          action2 = action2 + beta*(1.0_dp - real(plaquette(u,[x,y])))
       end do
    end do

    !action2 = beta*(L**2 - action2)
    
  end function action2
  
end module dynamics
