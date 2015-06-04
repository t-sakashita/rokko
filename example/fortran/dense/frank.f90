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
  type(rokko_serial_dense_solver) :: solver
  type(rokko_localized_matrix) :: mat, Z
  type(rokko_localized_vector) :: w
  character(len=100) :: solver_name, tmp_str
  integer arg_len, status

  integer :: provided, ierr, myrank, nprocs
  integer :: i

  if (command_argument_count().eq.2) then
     call get_command_argument(1, tmp_str, arg_len, status)
     solver_name = trim(tmp_str)
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     write(*,'(A)') "Error: frank solver_name dimension"
     stop
  endif

  write(*,*) "solver name = ", trim(solver_name)
  write(*,*) "matrix dimension = ", dim

  call rokko_serial_dense_solver_construct(solver, solver_name)

  call rokko_localized_matrix_construct(mat, dim, dim, rokko_matrix_col_major)
  call rokko_localized_matrix_construct(Z, dim, dim, rokko_matrix_col_major)
  call rokko_localized_vector_construct(w, dim)

  ! generate frank matrix
  call rokko_frank_matrix_generate_localized_matrix(mat)
  call rokko_localized_matrix_print(mat)

  call rokko_serial_dense_solver_diagonalize_localized_matrix(solver, mat, w, Z)

  write(*,'(A)') "Computed Eigenvalues = "
  do i = 1, dim
     write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
  enddo

  call rokko_localized_matrix_destruct(mat)
  call rokko_localized_matrix_destruct(Z)
  call rokko_localized_vector_destruct(w)
  call rokko_serial_dense_solver_destruct(solver)

end program frank_matrix
