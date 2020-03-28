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

  function rokko_get_string(s)
    use, intrinsic :: iso_c_binding
    interface
       pure function strlen(s) bind(c, name="strlen")
         use, intrinsic :: iso_c_binding
         type(c_ptr), intent(in), value :: s
         integer(c_size_t) :: strlen
       end function strlen
    end interface
    type(c_ptr), intent(in), value :: s
    character(kind=c_char, len=strlen(s)), pointer :: rokko_get_string

    call c_f_pointer(s, rokko_get_string)
  end function rokko_get_string

  subroutine string_free_c(str)
    type(c_ptr), intent(inout) :: str
    if (c_associated(str)) then
       call free_c(str)
       str = c_null_ptr
    end if
  end subroutine string_free_c
  
end module rokko_string
