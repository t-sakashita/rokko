!/*****************************************************************************
!*
!* Rokko: Integrated Interface for libraries of eigenvalue decomposition
!*
!* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!*                            Synge Todo <wistaria@comp-phys.org>,
!*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!*
!* Distributed under the Boost Software License, Version 1.0. (See accompanying
!* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!*
!*****************************************************************************/

program frank_matrix
  use rokko
  implicit none
  integer :: dim
  type(rokko_serial_dense_ev) :: solver
  type(rokko_localized_matrix) :: mat, Z
  type(rokko_localized_vector) :: w
  type(rokko_parameters) :: params, params_out
  character(len=:), allocatable :: library, routine
  character(len=20) :: solver_name, tmp_str
  integer arg_len, status

  integer :: provided, ierr, myrank, nprocs
  integer :: i, info

  if (command_argument_count() >= 1) then
     call get_command_argument(1, solver_name, arg_len, status)
  else
     solver_name = "lapack"
  endif
  call rokko_split_solver_name(solver_name, library, routine)

  if (command_argument_count() == 2) then  
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     dim = 10
  endif

  print *,"library = ", library
  print *,"routien = ", routine
  print *,"dimension = ", dim

  call rokko_serial_dense_ev_construct(solver, library)

  call rokko_localized_matrix_construct(mat, dim, dim, rokko_matrix_col_major)
  call rokko_localized_matrix_construct(Z, dim, dim, rokko_matrix_col_major)
  call rokko_localized_vector_construct(w, dim)

  ! generate frank matrix
  call rokko_frank_matrix_generate_localized_matrix(mat)
  call rokko_localized_matrix_print(mat)

  call rokko_parameters_construct(params)

  call rokko_parameters_set_string(params, "routine", "dsyevr")
  call rokko_parameters_set(params, "verbose", .true.)
  call rokko_parameters_set(params, "upper_index", 4)
  call rokko_parameters_set(params, "lower_index", 2)
  call rokko_parameters_set(params, "upper_value", 3.2d0)
  call rokko_parameters_set(params, "lower_value", 1.11d0)

!  call rokko_serial_dense_ev_diagonalize(solver, mat, w, Z, params, params_out)
!  call rokko_serial_dense_ev_diagonalize(solver, mat, w, Z, params)
!  call rokko_serial_dense_ev_diagonalize(solver, mat, w, Z)
  call rokko_serial_dense_ev_diagonalize(solver, mat, w, params, params_out)
!  call rokko_serial_dense_ev_diagonalize(solver, mat, w, params)
!  call rokko_serial_dense_ev_diagonalize(solver, mat, w)
  
  call rokko_parameters_get(params_out, "info", info)
!  call rokko_parameters_get_int(params_out, "info", info)

  print *, "info=", info

  write(*,'(A)') "Computed Eigenvalues = "
  do i = 1, dim
     write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
  enddo

  call rokko_localized_matrix_destruct(mat)
  call rokko_localized_matrix_destruct(Z)
  call rokko_localized_vector_destruct(w)
  call rokko_serial_dense_ev_destruct(solver)

end program frank_matrix
