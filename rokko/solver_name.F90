!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

#include <rokko/config.h>

module solver_name_utility
  use iso_c_binding
  implicit none

  interface

     ! interface for C function "void free(void *ptr)"
     subroutine free_c(ptr) bind(C,name="free")
       use iso_c_binding
       type(c_ptr), value, intent(in) :: ptr
     end subroutine free_c
     
     subroutine rokko_split_solver_name_f (str, library, routine) &
          bind(c,name='rokko_split_solver_name_f')
       use iso_c_binding
       implicit none
       character(c_char), intent(in) :: str(*)
       type(c_ptr), intent(out) :: library
       type(c_ptr), intent(out) :: routine
     end subroutine rokko_split_solver_name_f
     
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

  subroutine rokko_split_solver_name(str, library, routine)
    character(*), intent(in) :: str
    character(len=:), allocatable, intent(out) :: library
    character(len=:), allocatable, intent(out) :: routine
    type(c_ptr) :: library_ptr, routine_ptr
    character*255 :: tmp
    call rokko_split_solver_name_f (str//c_null_char, library_ptr, routine_ptr)
    call rokko_get_string (library_ptr, library)
    call rokko_get_string (routine_ptr, routine)
  end subroutine rokko_split_solver_name

end module solver_name_utility


