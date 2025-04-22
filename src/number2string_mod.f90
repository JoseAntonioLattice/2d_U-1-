module number2string_mod

    use iso_fortran_env, only : dp => real64, i4 => int32

    implicit none

    private

    public :: int2str, real2str

contains

    function int2str(k)
      !"Convert an integer to string."
      integer(i4), intent(in) :: k
      character(:), allocatable :: int2str
      character(20) :: var
      write (var, *) k
      int2str = trim(adjustl(var))
    end function int2str

    character(6) function real2str(real_)
      real(dp), intent(in) :: real_
      character(6) :: string_
      write(string_,"(f6.4)") real_
      real2str = string_
    end function real2str

end module number2string_mod
