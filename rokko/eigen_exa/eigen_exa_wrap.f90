!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine eigen_init_wrap(comm, order)
  use eigen_libs, only : eigen_init
  implicit none
  integer, intent(in), optional :: comm
  character*(*), intent(in), optional :: order
  call eigen_init( comm, order )
  return
end subroutine eigen_init_wrap

subroutine eigen_free_wrap(flag)
  use eigen_libs, only : eigen_free
  implicit none
  integer, intent(in), optional ::  flag
  call eigen_free(flag)
  return
end subroutine eigen_free_wrap
