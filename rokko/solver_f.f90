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

subroutine rokko_solver_construct(solver, solver_name)
  use iso_c_binding
  use rokko, only: rokko_solver
  implicit none
  interface
     subroutine rokko_solver_construct_f(solver, solver_name) bind(c)
       use iso_c_binding
       use rokko, only: rokko_solver
       implicit none
       type(rokko_solver), intent(out) :: solver
       character(kind=c_char), intent(in) :: solver_name(*)
     end subroutine rokko_solver_construct_f
  end interface
  type(rokko_solver), intent(inout) :: solver
  character(*), intent(in) :: solver_name
  call rokko_solver_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
end subroutine rokko_solver_construct
