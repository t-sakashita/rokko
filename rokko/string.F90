!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

module rokko_string
  use iso_c_binding
  implicit none

  type string
     character(len=:), allocatable :: str
  end type string

  type array_strings
     type(string), allocatable :: string(:)
     integer :: size
  end type array_strings
  
  interface
     ! interface for C function "void free(void *ptr)"
     subroutine free_c(ptr) bind(C,name="free")
       use iso_c_binding
       type(c_ptr), value, intent(in) :: ptr
     end subroutine free_c

     type(c_ptr) function rokko_string_i_c (ptr, i) &
          bind(c,name='rokko_string_i')
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ptr
       integer(c_int), value, intent(in) :: i
     end function rokko_string_i_c
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

  subroutine rokko_get_string_fixedsize (str_ptr, val)
    use iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: str_ptr
    character(len=*), intent(out) :: val
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
  end subroutine rokko_get_string_fixedsize
  
  subroutine string_free_c(str)
    use iso_c_binding
    type(c_ptr), intent(inout) :: str
    if (c_associated(str)) then
       call free_c(str)
       str = c_null_ptr
    end if
  end subroutine string_free_c
  
end module rokko_string
