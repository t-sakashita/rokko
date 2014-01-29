!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>,
!                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine rokko_grid_construct(grid, comm, grid_major)
  use iso_c_binding
  use rokko_f, only: rokko_grid, rokko_grid_col_major
  implicit none
  interface
     subroutine rokko_grid_construct_f(grid, comm, grid_major) bind(c)
       use iso_c_binding
       use rokko_f, only: rokko_grid
       implicit none
       type(rokko_grid), intent(out) :: grid
       integer(c_int), value, intent(in) :: comm
       integer(c_int), value, intent(in) :: grid_major
     end subroutine rokko_grid_construct_f
  end interface
  type(rokko_grid), intent(out) :: grid
  integer, intent(in) :: comm
  integer, intent(in) :: grid_major
  integer(c_int) :: comm_out
  integer(c_int) :: grid_major_out
  if (present(grid_major)) then
     grid_major_out = grid_major
  else
     grid_major_out = rokko_grid_col_major
  endif
  call rokko_grid_construct_f(grid, comm_out, grid_major_out)
end subroutine rokko_grid_construct
