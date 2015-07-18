module rokko_string
  use iso_c_binding
  implicit none

  interface

     ! interface for C function "void free(void *ptr)"
     subroutine free_c(ptr) bind(C,name="free")
       use iso_c_binding
       type(c_ptr), value, intent(in) :: ptr
     end subroutine free_c
     
  end interface

contains

  subroutine rokko_get_string (str_ptr, val)
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: str_ptr
    character(len=:), allocatable, intent(out) :: val
    character, pointer, dimension(:) :: tmp_array
    character*255 :: tmp
    integer :: i, n
    call c_f_pointer(str_ptr, tmp_array, [ 255 ] )
    do i=1, 255
       if (tmp_array(i) == char(0)) exit
       tmp(i:i) = tmp_array(i)
    enddo
    n = i - 1
    call free_c(str_ptr)
    val = trim(tmp(1:n))  ! automatically allocating suitable size
  end subroutine rokko_get_string

  subroutine string_free_c(str)
    use iso_c_binding
    type(c_ptr), intent(inout) :: str
    if (c_associated(str)) then
       call free_c(str)
       str = c_null_ptr
    end if
  end subroutine string_free_c

end module rokko_string
