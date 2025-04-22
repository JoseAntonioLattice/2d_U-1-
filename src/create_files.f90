module create_files
  use check_files_directories_mod
  use number2string_mod
  implicit none

contains

  subroutine create_measurements_file(L,DATA_FILE)
    integer, intent(in) :: L
    !real(8), intent(in) :: beta
    !character(*), intent(in) :: algorithm
    logical :: file_exists, condition
    integer :: i

    character(100) :: directory
    character(:), allocatable, intent(out) :: data_file

    
    directory = "data"
    call check_directory(trim(directory))
    
    directory = trim(directory)//"/L="//trim(int2str(L))
    call check_directory(trim(directory))

    directory = trim(directory)//"/measurements_L="//trim(int2str(L))//"_"

    i = 1
    do 
       data_file = trim(directory)//trim(int2str(i))//".dat"
       call check_file(trim(data_file), file_exists)
       condition = file_exists
       if (condition .eqv. .false.) exit
       i = i + 1
    end do
    !open(unit = 100, file = data_file)
    !write(100, nml = input_parameters )
    
  end subroutine create_measurements_file

end module create_files
