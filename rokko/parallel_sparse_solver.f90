!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine rokko_parallel_sparse_solver_construct(solver, solver_name)
  use iso_c_binding
  use rokko, only: rokko_parallel_sparse_solver
  implicit none
  interface
     subroutine rokko_parallel_sparse_solver_construct_f(solver, solver_name) bind(c)
       use rokko
       implicit none
       type(rokko_parallel_sparse_solver), intent(out) :: solver
       character(kind=c_char), intent(in) :: solver_name(*)
     end subroutine rokko_parallel_sparse_solver_construct_f
  end interface
  type(rokko_parallel_sparse_solver), intent(inout) :: solver
  character(*), intent(in) :: solver_name
  call rokko_parallel_sparse_solver_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
end subroutine rokko_parallel_sparse_solver_construct
