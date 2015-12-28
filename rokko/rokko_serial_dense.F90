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

#include <rokko/config.h>

module rokko_serial_dense
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_grid_col_major = 1, rokko_grid_row_major = 2
     enumerator :: rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
  end enum

  !
  ! classes
  !

  type, bind(c) :: rokko_serial_dense_ev
     type(c_ptr) ptr
  end type rokko_serial_dense_ev

  type, bind(c) :: rokko_localized_vector
     type(c_ptr) ptr
  end type rokko_localized_vector

  type, bind(c) :: rokko_localized_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_localized_matrix
  
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
     function rokko_localized_vector_get_f(vec, i) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       real(c_double) :: rokko_localized_vector_get_f
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
  ! rokko_serial_dense_ev
  !

  interface
     subroutine rokko_serial_dense_ev_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_serial_dense_ev
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
     end subroutine rokko_serial_dense_ev_destruct
  end interface

  interface rokko_serial_dense_ev_diagonalize

     subroutine rokko_serial_dense_ev_diagonalize(solver, mat, &
          eigvals, eigvecs, params, params_out) bind(c,name="rokko_serial_dense_ev_diagonalize_f")
       use iso_c_binding
       use parameters
       import rokko_serial_dense_ev, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_localized_matrix), intent(inout) :: eigvecs
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_serial_dense_ev_diagonalize
     
     subroutine rokko_serial_dense_ev_diagonalize_no_params_out(solver, mat, &
          eigvals, eigvecs, params) bind(c,name="rokko_serial_dense_ev_diagonalize_no_params_out_f")
       use iso_c_binding
       use parameters
       import rokko_serial_dense_ev, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_localized_matrix), intent(inout) :: eigvecs
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_serial_dense_ev_diagonalize_no_params_out
          
     subroutine rokko_serial_dense_ev_diagonalize_no_params_inout(solver, mat, &
          eigvals, eigvecs) bind(c,name='rokko_serial_dense_ev_diagonalize_no_params_inout_f')
       use iso_c_binding
       use parameters
       import rokko_serial_dense_ev, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_localized_matrix), intent(inout) :: eigvecs
     end subroutine rokko_serial_dense_ev_diagonalize_no_params_inout

     ! Only eigenvalues
     subroutine rokko_serial_dense_ev_diagonalize_eigvals(solver, mat, eigvals, params, params_out) &
          bind(c,name="rokko_serial_dense_ev_diagonalize_eigvals_f")
       use iso_c_binding
       use parameters
       import rokko_serial_dense_ev, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_serial_dense_ev_diagonalize_eigvals
     
     subroutine rokko_serial_dense_ev_diagonalize_eigvals_no_params_out(solver, mat, eigvals, params) &
          bind(c,name="rokko_serial_dense_ev_diagonalize_eigvals_no_params_out_f")
       use iso_c_binding
       use parameters
       import rokko_serial_dense_ev, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_serial_dense_ev_diagonalize_eigvals_no_params_out
     
     subroutine rokko_serial_dense_ev_diagonalize_eigvals_no_params_inout(solver, mat, eigvals) &
          bind(c,name="rokko_serial_dense_ev_diagonalize_eigvals_no_params_inout_f")
       use iso_c_binding
       use parameters
       import rokko_serial_dense_ev, rokko_localized_matrix, rokko_localized_vector
       implicit none
       type(rokko_serial_dense_ev), intent(inout) :: solver
       type(rokko_localized_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
     end subroutine rokko_serial_dense_ev_diagonalize_eigvals_no_params_inout     
     
  end interface rokko_serial_dense_ev_diagonalize

  !
  ! rokko_frank_matrix for serial solvers
  !

  interface
     subroutine rokko_frank_matrix_generate_localized_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), intent(inout) :: matrix
     end subroutine rokko_frank_matrix_generate_localized_matrix
  end interface

contains

  subroutine rokko_serial_dense_ev_construct(solver, solver_name)
    use iso_c_binding
    implicit none
    interface
       subroutine rokko_serial_dense_ev_construct_f(solver, solver_name) bind(c)
         use iso_c_binding
         import rokko_serial_dense_ev
         implicit none
         type(rokko_serial_dense_ev), intent(out) :: solver
         character(kind=c_char), intent(in) :: solver_name(*)
       end subroutine rokko_serial_dense_ev_construct_f
    end interface
    type(rokko_serial_dense_ev), intent(inout) :: solver
    character(*), intent(in) :: solver_name
    call rokko_serial_dense_ev_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
  end subroutine rokko_serial_dense_ev_construct

end module rokko_serial_dense

