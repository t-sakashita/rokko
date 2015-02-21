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

module rokko
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_grid_col_major = 1, rokko_grid_row_major = 2
     enumerator :: rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
  end enum

  !
  ! classes
  !
  
  type, bind(c) :: rokko_grid
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_grid

  type, bind(c) :: rokko_serial_dense_solver
     type(c_ptr) ptr
  end type rokko_serial_dense_solver

  type, bind(c) :: rokko_parallel_dense_solver
     type(c_ptr) ptr
  end type rokko_parallel_dense_solver
  
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

  type, bind(c) :: rokko_distributed_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_distributed_matrix

  type, bind(c) :: rokko_distributed_crs_matrix
     type(c_ptr) ptr
  end type rokko_distributed_crs_matrix

  !
  ! rokko_grid
  !
  
  interface
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
  end interface
  
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

  !
  ! rokko_serial_dense_solver
  !

  interface
     subroutine rokko_serial_dense_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_serial_dense_solver
       implicit none
       type(rokko_serial_dense_solver), intent(inout) :: solver
     end subroutine rokko_serial_dense_solver_destruct
     
     subroutine rokko_serial_dense_solver_diagonalize_localized_matrix(solver, mat, eigvals, &
          eigvecs) bind(c)
       use iso_c_binding
       import rokko_serial_dense_solver, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_solver), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_localized_matrix), intent(inout) :: eigvecs
     end subroutine rokko_serial_dense_solver_diagonalize_localized_matrix
  end interface
  
  !
  ! rokko_parallel_dense_solver
  !
     
  interface
     subroutine rokko_parallel_dense_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_dense_solver
       implicit none
       type(rokko_parallel_dense_solver), intent(inout) :: solver
     end subroutine rokko_parallel_dense_solver_destruct
     
     subroutine rokko_parallel_dense_solver_diagonalize_distributed_matrix(solver, mat, eigvals, &
          eigvecs) bind(c)
       use iso_c_binding
       import rokko_parallel_dense_solver, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_solver), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_matrix), intent(inout) :: eigvecs
     end subroutine rokko_parallel_dense_solver_diagonalize_distributed_matrix
  end interface

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
  ! rokko_localized_vector
  !

  interface
     subroutine rokko_localized_vector_construct(vec, dim1) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), intent(out) :: vec
       integer(c_int), value, intent(in) :: dim1
     end subroutine rokko_localized_vector_construct

     subroutine rokko_localized_vector_destruct(vec) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), intent(inout) :: vec
     end subroutine rokko_localized_vector_destruct
  end interface

  interface rokko_localized_vector_get
     real(c_double) function rokko_localized_vector_get_f(vec, i) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
     end function rokko_localized_vector_get_f
  end interface rokko_localized_vector_get

  !
  ! rokko_localized_matrix
  !

  interface
     subroutine rokko_localized_matrix_construct(matrix, dim1, dim2, matrix_major) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       integer(c_int), value, intent(in) :: matrix_major
     end subroutine rokko_localized_matrix_construct
     
     subroutine rokko_localized_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), intent(inout) :: matrix
     end subroutine rokko_localized_matrix_destruct

     subroutine rokko_localized_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), value, intent(in) :: matrix
     end subroutine rokko_localized_matrix_print
  end interface
  
  !
  ! rokko_distributed_matrix
  !

  interface
     subroutine rokko_distributed_matrix_construct(matrix, dim1, dim2, grid, solver, matrix_major) &
          bind(c)
       use iso_c_binding
       import rokko_grid, rokko_parallel_dense_solver, rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       type(rokko_grid), value, intent(in) :: grid
       type(rokko_parallel_dense_solver), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: matrix_major
     end subroutine rokko_distributed_matrix_construct
     
     subroutine rokko_distributed_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_matrix_destruct
     
     subroutine rokko_distributed_matrix_generate_array(matrix, array)
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       real(8), intent(in) :: array(:,:)
     end subroutine rokko_distributed_matrix_generate_array

     subroutine rokko_distributed_matrix_generate_function(matrix, func) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(c_funptr), intent(in), value :: func
     end subroutine rokko_distributed_matrix_generate_function
     
     subroutine rokko_distributed_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_matrix_print
     
     subroutine rokko_distributed_matrix_set_local(matrix, local_i, local_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: local_i, local_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_local
     
     real(c_double) function rokko_distributed_matrix_get_local(matrix, local_i,local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i,local_j
     end function rokko_distributed_matrix_get_local
     
     subroutine rokko_distributed_matrix_set_global(matrix, global_i, global_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: global_i, global_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_global
     
     real(c_double) function rokko_distributed_matrix_get_global(matrix, global_i, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in):: global_i, global_j
     end function rokko_distributed_matrix_get_global
     
     integer(c_int) function rokko_distributed_matrix_get_m_local(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_local

     integer(c_int) function rokko_distributed_matrix_get_n_local(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_local

     integer(c_int) function rokko_distributed_matrix_get_m_global(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_global

     integer(c_int) function rokko_distributed_matrix_get_n_global(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_global

     integer(c_int) function rokko_distributed_matrix_get_nprocs(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_nprocs

     integer(c_int) function rokko_distributed_matrix_get_myrank(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_myrank
 
     integer(c_int) function rokko_distributed_matrix_translate_l2g_row(matrix, local_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value,intent(in) :: local_i
     end function rokko_distributed_matrix_translate_l2g_row

     integer(c_int) function rokko_distributed_matrix_translate_l2g_col(matrix, local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::local_j
     end function rokko_distributed_matrix_translate_l2g_col

     integer(c_int) function rokko_distributed_matrix_translate_g2l_row(matrix, global_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::global_i
     end function rokko_distributed_matrix_translate_g2l_row

     integer(c_int) function rokko_distributed_matrix_translate_g2l_col(matrix, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_j
     end function rokko_distributed_matrix_translate_g2l_col
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
     
     integer(c_int) function rokko_distributed_crs_matrix_start_row(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_start_row
     
     integer(c_int) function rokko_distributed_crs_matrix_end_row(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
     end function rokko_distributed_crs_matrix_end_row
     
     integer(c_int) function rokko_distributed_crs_matrix_num_local_rows(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
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
  ! collective operations
  !

  interface
     integer(c_int) function rokko_gather(matrix, array, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value ::root
     end function rokko_gather

     subroutine rokko_all_gather(matrix, array)
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       real(8), intent(in), target :: array(:,:)
     end subroutine rokko_all_gather

     integer(c_int) function rokko_scatter(array, matrix, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value :: root
     end function rokko_scatter
  end interface

  !
  ! rokko_frank_matrix
  !

  interface
     subroutine rokko_frank_matrix_generate_localized_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), intent(inout) :: matrix
     end subroutine rokko_frank_matrix_generate_localized_matrix

     subroutine rokko_frank_matrix_generate_distributed_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_frank_matrix_generate_distributed_matrix
  end interface
  
contains

  subroutine rokko_serial_dense_solver_construct(solver, solver_name)
    use iso_c_binding
    implicit none
    interface
       subroutine rokko_serial_dense_solver_construct_f(solver, solver_name) bind(c)
         use iso_c_binding
         import rokko_serial_dense_solver
         implicit none
         type(rokko_serial_dense_solver), intent(out) :: solver
         character(kind=c_char), intent(in) :: solver_name(*)
       end subroutine rokko_serial_dense_solver_construct_f
    end interface
    type(rokko_serial_dense_solver), intent(inout) :: solver
    character(*), intent(in) :: solver_name
    call rokko_serial_dense_solver_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
  end subroutine rokko_serial_dense_solver_construct

end module rokko
