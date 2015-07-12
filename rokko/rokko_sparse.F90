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

module rokko_sparse
  use iso_c_binding
  implicit none

  !
  ! classes
  !
  
  type, bind(c) :: rokko_parallel_sparse_solver
     type(c_ptr) ptr
  end type rokko_parallel_sparse_solver

  type, bind(c) :: rokko_distributed_crs_matrix
     type(c_ptr) ptr
  end type rokko_distributed_crs_matrix

  type, bind(c) :: rokko_distributed_mfree
     type(c_ptr) ptr
  end type rokko_distributed_mfree

  !
  ! rokko_parallel_sparse_solver
  !

  interface
     subroutine rokko_parallel_sparse_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
     end subroutine rokko_parallel_sparse_solver_destruct

     subroutine rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(solver, mat, &
          num_evals, block_size, max_iters, tol) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       integer(c_int), value, intent(in) :: num_evals, block_size, max_iters
       real(8), value, intent(in) :: tol
     end subroutine rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix

     subroutine rokko_parallel_sparse_solver_diagonalize_distributed_mfree(solver, mat, &
          num_evals, block_size, max_iters, tol) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_mfree
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       type(rokko_distributed_mfree), intent(inout) :: mat
       integer(c_int), value, intent(in) :: num_evals, block_size, max_iters
       real(8), value, intent(in) :: tol
     end subroutine rokko_parallel_sparse_solver_diagonalize_distributed_mfree

     function rokko_parallel_sparse_solver_num_conv(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_parallel_sparse_solver_num_conv
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
     end function rokko_parallel_sparse_solver_num_conv

     function rokko_parallel_sparse_solver_eigenvalue(solver, i) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       real(c_double) :: rokko_parallel_sparse_solver_eigenvalue
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       integer(c_int), value, intent(in) :: i
     end function rokko_parallel_sparse_solver_eigenvalue

     subroutine rokko_parallel_sparse_solver_eigenvector(solver, i, eig_vec) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_solver, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_solver), intent(inout) :: solver
       integer(c_int), value, intent(in) :: i
       real(c_double), dimension(*), intent(inout) :: eig_vec
     end subroutine rokko_parallel_sparse_solver_eigenvector
  end interface

  !
  ! rokko_distributed_crs_matrix
  !

  interface
     subroutine rokko_distributed_crs_matrix_construct(matrix, dim1, dim2, solver) bind(c)
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
       integer(c_int), dimension(col_size), intent(in) :: cols
       real(c_double), dimension(col_size), intent(in) :: values
     end subroutine rokko_distributed_crs_matrix_insert
     
     subroutine rokko_distributed_crs_matrix_complete(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_crs_matrix_complete
     
     function rokko_distributed_crs_matrix_start_row(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_start_row
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_start_row
     
     function rokko_distributed_crs_matrix_end_row(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_end_row
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_end_row
     
     function rokko_distributed_crs_matrix_num_local_rows(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_num_local_rows
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_num_local_rows
     
     subroutine rokko_distributed_crs_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end subroutine rokko_distributed_crs_matrix_print
  end interface

  !
  ! rokko_distributed_mfree
  !

  interface
     subroutine rokko_distributed_mfree_f_construct(matrix, func, dim, num_local_rows) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(out) :: matrix
       type(c_funptr), intent(in), value :: func
       integer(c_int), value, intent(in) :: dim, num_local_rows
     end subroutine rokko_distributed_mfree_f_construct
    
     integer(c_int) function rokko_distributed_mfree_num_local_rows(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(in) :: matrix
     end function rokko_distributed_mfree_num_local_rows
     
     integer(c_int) function rokko_distributed_mfree_dim(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(in) :: matrix
     end function rokko_distributed_mfree_dim
     
     subroutine rokko_distributed_mfree_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(inout) :: matrix
     end subroutine rokko_distributed_mfree_destruct     
  end interface

contains
  SUBROUTINE rokko_distributed_mfree_construct(mat, multiply_in, dim, num_local_rows)
    USE, INTRINSIC :: ISO_C_BINDING
    type(rokko_distributed_mfree), intent(inout) :: mat
    INTEGER(C_INT), INTENT(IN) :: dim, num_local_rows
    TYPE(C_FUNPTR) :: cproc
    INTERFACE
       SUBROUTINE multiply_in (n, x, y) BIND(C)
         USE, INTRINSIC :: ISO_C_BINDING
         INTEGER(C_INT), INTENT(IN), VALUE :: n
         REAL(C_DOUBLE), INTENT(IN) :: x(n)
         REAL(C_DOUBLE), INTENT(OUT) :: y(n)
       END SUBROUTINE multiply_in
    END INTERFACE
    ! Get C procedure pointer.
    cproc = C_FUNLOC(multiply_in)
    ! call wrapper written in C.
    call rokko_distributed_mfree_f_construct(mat, cproc, dim, num_local_rows)
  END SUBROUTINE rokko_distributed_mfree_construct

end module rokko_sparse

