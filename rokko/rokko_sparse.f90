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

  type, bind(c) :: rokko_parallel_sparse_solver
     type(c_ptr) ptr
  end type rokko_parallel_sparse_solver
  
  type, bind(c) :: rokko_distributed_crs_matrix
     type(c_ptr) ptr
  end type rokko_distributed_crs_matrix

  interface

     ! rokko_parallel_sparse_solver
     
     subroutine rokko_parallel_sparse_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
     end subroutine rokko_parallel_sparse_solver_destruct

     subroutine rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(solver, mat, &
          & num_evals, block_size, max_iters, tol) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       integer(c_int), value, intent(in) :: num_evals, block_size, max_iters
       real(8), value, intent(in) :: tol
     end subroutine rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix

     integer(c_int) function rokko_parallel_sparse_solver_num_conv(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
     end function rokko_parallel_sparse_solver_num_conv

     real(c_double) function rokko_parallel_sparse_solver_eigenvalue(solver, i) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       integer(c_int), value, intent(in) :: i
     end function rokko_parallel_sparse_solver_eigenvalue

     real(c_double) function rokko_parallel_sparse_solver_eigenvector(solver, i) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       integer(c_int), value, intent(in) :: i
     end function rokko_parallel_sparse_solver_eigenvector

     
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

     subroutine rokko_distributed_crs_matrix_insert(matrix, row, col_size, cols, values) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(inout) :: matrix
       integer(c_int), value, intent(in) :: row, col_size
       integer(c_int), dimension(col_size) :: cols
       real(c_double), dimension(col_size) :: values
     end subroutine rokko_distributed_crs_matrix_insert

     subroutine rokko_distributed_crs_matrix_complete(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_crs_matrix_complete

     integer(c_int) function rokko_distributed_crs_matrix_start_row(matrix) &
          bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_start_row

     integer(c_int) function rokko_distributed_crs_matrix_end_row(matrix) &
          bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_end_row

     integer(c_int) function rokko_distributed_crs_matrix_num_local_rows(matrix) &
          bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_num_local_rows

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
