!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>,
!                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>,
!                            Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module rokko_sparse
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_grid_col_major = 1, rokko_grid_row_major = 2
     enumerator :: rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
  end enum

  type, bind(c) :: rokko_grid
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_grid

  type, bind(c) :: rokko_parallel_sparse_solver
     type(c_ptr) ptr
  end type rokko_parallel_sparse_solver
  
  type, bind(c) :: rokko_localized_vector
     type(c_ptr) ptr
  end type rokko_localized_vector

  type, bind(c) :: rokko_localized_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_localized_matrix

  type, bind(c) :: rokko_distributed_crs_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_distributed_crs_matrix

  interface rokko_grid_construct
     subroutine rokko_grid_construct_f(grid, comm, grid_major) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       type(rokko_grid), intent(out) :: grid
       integer(c_int), value, intent(in) :: comm
       integer(c_int), value, intent(in) :: grid_major
     end subroutine rokko_grid_construct_f
  end interface rokko_grid_construct
  
  interface

     ! rokko_grid
     
     subroutine rokko_grid_destruct(grid) bind(c)
       import rokko_grid
       implicit none
       type(rokko_grid), intent(inout) :: grid
     end subroutine rokko_grid_destruct

     integer(c_int) function rokko_grid_get_myrank(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_myrank

     integer(c_int) function rokko_grid_get_nprocs(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_nprocs

     ! rokko_parallel_sparse_solver
     
     subroutine rokko_parallel_sparse_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
     end subroutine rokko_parallel_sparse_solver_destruct
     
     subroutine rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(solver, mat, eigvals, eigvecs) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_crs_matrix), intent(inout) :: eigvecs
     end subroutine rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix
     
     ! rokko_distributed_crs_matrix

     subroutine rokko_distributed_crs_matrix_construct(matrix, dim1, dim2, solver) &
          bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       type(rokko_parallel_sparse_solver), value, intent(in) :: solver
     end subroutine rokko_distributed_crs_matrix_construct
     
     subroutine rokko_distributed_crs_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_crs_matrix_destruct

  end interface
contains

subroutine rokko_parallel_sparse_solver_construct(solver, solver_name)
  use iso_c_binding
!  import rokko_parallel_sparse_solver
  implicit none
  interface
     subroutine rokko_parallel_sparse_solver_construct_f(solver, solver_name) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver
       implicit none
       type(rokko_parallel_sparse_solver), intent(out) :: solver
       character(kind=c_char), intent(in) :: solver_name(*)
     end subroutine rokko_parallel_sparse_solver_construct_f
  end interface
  type(rokko_parallel_sparse_solver), intent(inout) :: solver
  character(*), intent(in) :: solver_name
  call rokko_parallel_sparse_solver_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
end subroutine rokko_parallel_sparse_solver_construct

end module rokko_sparse
