module pbc

  use iso_fortran_env, only : i4 => int32
  
  implicit none

  integer, allocatable, dimension(:) :: ip1, im1, ip2, im2

contains

  subroutine set_pbc(L)
    
    integer(i4), intent(in) :: L(2)

    allocate(ip1(L(1)), im1(L(1)))
    allocate(ip2(L(2)), im2(L(2)))
    call initialize(L)
    
  end subroutine set_pbc

  subroutine initialize(L)

    integer(i4), intent(in) :: L(2)
    integer(i4) :: i

    do i = 1, L(1)
       ip1(i) = i + 1
       im1(i) = i - 1
    end do
    ip1(L(1)) = 1
    im1(1) = L(1)

    
    do i = 1, L(2)
       ip2(i) = i + 1
       im2(i) = i - 1
    end do
    ip2(L(2)) = 1
    im2(1) = L(2)
    
  end subroutine initialize

  function ipf(vector, mu)

    integer(i4), dimension(:), intent(in) :: vector
    integer(i4) :: mu

    integer(i4), dimension(size(vector)) :: ipf
    
    ipf = vector

    select case(mu)
    case(1)
       ipf(mu) = ip1(vector(mu))
    case(2)
       ipf(mu) = ip2(vector(mu))
    end select
  end function ipf

  function imf(vector, mu)

    integer(i4), dimension(:), intent(in) :: vector
    integer(i4) :: mu

    integer(i4), dimension(size(vector)) :: imf
    
    imf = vector

    select case(mu)
    case(1)
       imf(mu) = im1(vector(mu))
    case(2)
       imf(mu) = im2(vector(mu))
    end select
    
  end function imf
  
end module pbc
