!*****************************************************************************
!
! rokko: integrated interface for libraries of eigenvalue decomposition
!
! copyright (c) 2012-2015 by rokko developers https://github.com/t-sakashita/rokko
!
! distributed under the boost software license, version 1.0. (see accompanying
! file license_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

module rokko_sparse
  use iso_c_binding
  implicit none

  !
  ! classes
  !
  
  type, bind(c) :: rokko_parallel_sparse_ev
     type(c_ptr) ptr
  end type rokko_parallel_sparse_ev

  type, bind(c) :: rokko_distributed_crs_matrix
     type(c_ptr) ptr
  end type rokko_distributed_crs_matrix

  type, bind(c) :: rokko_distributed_mfree
     type(c_ptr) ptr
  end type rokko_distributed_mfree

  !
  ! rokko_parallel_sparse_ev
  !

  interface
     subroutine rokko_parallel_sparse_ev_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
     end subroutine rokko_parallel_sparse_ev_destruct

     function rokko_parallel_sparse_ev_num_conv(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_parallel_sparse_ev_num_conv
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
     end function rokko_parallel_sparse_ev_num_conv

     function rokko_parallel_sparse_ev_eigenvalue(solver, i) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       real(c_double) :: rokko_parallel_sparse_ev_eigenvalue
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: i
     end function rokko_parallel_sparse_ev_eigenvalue

     subroutine rokko_parallel_sparse_ev_eigenvector(solver, i, eig_vec) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: i
       real(c_double), dimension(*), intent(inout) :: eig_vec
     end subroutine rokko_parallel_sparse_ev_eigenvector
  end interface

  interface
     type(c_ptr) function rokko_parallel_sparse_ev_default_solver_c() &
          bind(c,name='rokko_parallel_sparse_ev_default_solver')
       use iso_c_binding
       implicit none
     end function rokko_parallel_sparse_ev_default_solver_c
  end interface
  
  interface rokko_parallel_sparse_ev_diagonalize
      
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(solver, mat, params, params_out) &
          bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix
     
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_noreturn(solver, &
          mat, params) bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_noreturn_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_noreturn
     
     
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree(solver, mat, params, params_out) &
          bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_mfree_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_mfree
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_mfree), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree
     
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn(solver, mat, params) &
          bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_mfree
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_mfree), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn

  end interface rokko_parallel_sparse_ev_diagonalize

  !
  ! rokko_distributed_crs_matrix
  !

  interface
     subroutine rokko_distributed_crs_matrix_construct(matrix, dim1, dim2, solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
     end subroutine rokko_distributed_crs_matrix_construct
     
     subroutine rokko_distributed_crs_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_crs_matrix_destruct
     
     subroutine rokko_distributed_crs_matrix_insert_c(matrix, row, col_size, cols, values) &
          bind(c,name='rokko_distributed_crs_matrix_insert')
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: row, col_size
       integer(c_int), dimension(col_size), intent(in) :: cols
       real(c_double), dimension(col_size), intent(in) :: values
     end subroutine rokko_distributed_crs_matrix_insert_c
     
     subroutine rokko_distributed_crs_matrix_complete(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_crs_matrix_complete
     
     function rokko_distributed_crs_matrix_start_row_c(matrix) &
          & bind(c,name='rokko_distributed_crs_matrix_start_row')
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_start_row_c
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end function rokko_distributed_crs_matrix_start_row_c
     
     function rokko_distributed_crs_matrix_end_row_c(matrix) &
       bind(c,name='rokko_distributed_crs_matrix_end_row')
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_end_row_c
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end function rokko_distributed_crs_matrix_end_row_c
     
     function rokko_distributed_crs_matrix_num_local_rows(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_num_local_rows
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end function rokko_distributed_crs_matrix_num_local_rows
     
     subroutine rokko_distributed_crs_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
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
       type(rokko_distributed_mfree), value, intent(in) :: matrix
     end function rokko_distributed_mfree_num_local_rows
     
     integer(c_int) function rokko_distributed_mfree_dim(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), value, intent(in) :: matrix
     end function rokko_distributed_mfree_dim
     
     subroutine rokko_distributed_mfree_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(inout) :: matrix
     end subroutine rokko_distributed_mfree_destruct     
  end interface

contains

  subroutine rokko_parallel_sparse_ev_default_solver(name)
    use rokko_string
    character(len=*), intent(out) :: name
    type(c_ptr) :: name_ptr
    name_ptr = rokko_parallel_sparse_ev_default_solver_c ()
    call rokko_get_string_fixedsize (name_ptr, name)
  end subroutine rokko_parallel_sparse_ev_default_solver
  
  subroutine rokko_distributed_mfree_construct(mat, multiply_in, dim, num_local_rows)
    use, intrinsic :: iso_c_binding
    type(rokko_distributed_mfree), intent(inout) :: mat
    integer(c_int), intent(in) :: dim, num_local_rows
    type(c_funptr) :: cproc
    interface
       subroutine multiply_in (n, x, y) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int), intent(in), value :: n
         real(c_double), intent(in) :: x(n)
         real(c_double), intent(out) :: y(n)
       end subroutine multiply_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(multiply_in)
    ! call wrapper written in c.
    call rokko_distributed_mfree_f_construct(mat, cproc, dim, num_local_rows)
  end subroutine rokko_distributed_mfree_construct

end module rokko_sparse

