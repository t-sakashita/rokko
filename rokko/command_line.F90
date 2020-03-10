!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module command_line_mod
  implicit none

contains

  subroutine get_command_argument_deferred(i, str)
    implicit none
    integer, intent(in) :: i
    character(len=:), allocatable, intent(out) :: str
    integer :: arg_len

    call get_command_argument(i, length=arg_len)
    allocate(character(arg_len) :: str)
    call get_command_argument(i, value=str)
  end subroutine get_command_argument_deferred

end module command_line_mod
