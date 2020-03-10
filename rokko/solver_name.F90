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
  use rokko_string
  implicit none

  interface
     
     subroutine rokko_split_solver_name_f (str, library, routine) &
          & bind(c,name='rokko_split_solver_name')
       use iso_c_binding
       implicit none
       character(c_char), intent(in) :: str(*)
       type(c_ptr), intent(out) :: library
       type(c_ptr), intent(out) :: routine
     end subroutine rokko_split_solver_name_f
     
  end interface

contains

  subroutine rokko_split_solver_name(str, library, routine)
    character(len=:), allocatable, intent(in) :: str
    character(len=:), allocatable, intent(out) :: library
    character(len=:), allocatable, intent(out) :: routine
    type(c_ptr) :: library_ptr, routine_ptr
    call rokko_split_solver_name_f (str//c_null_char, library_ptr, routine_ptr)
    call rokko_get_string (library_ptr, library)
    call rokko_get_string (routine_ptr, routine)
  end subroutine rokko_split_solver_name

  subroutine rokko_split_solver_name_fixedsize(str, library, routine)
    character(*), intent(in) :: str
    character(len=*), intent(out) :: library
    character(len=*), intent(out) :: routine
    type(c_ptr) :: library_ptr, routine_ptr
    call rokko_split_solver_name_f (str//c_null_char, library_ptr, routine_ptr)
    call rokko_get_string_fixedsize (library_ptr, library)
    call rokko_get_string_fixedsize (routine_ptr, routine)
  end subroutine rokko_split_solver_name_fixedsize

end module solver_name_utility


