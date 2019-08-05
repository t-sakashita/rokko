!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program frank_matrix
  use rokko
  implicit none
  integer :: dim
  type(rokko_serial_dense_ev) :: solver
  type(rokko_localized_matrix) :: mat, Z
  type(rokko_localized_vector) :: w
  type(rokko_parameters) :: params, params_out
  character(len=20) :: library, routine
  character(len=20) :: solver_name, tmp_str
  integer arg_len, status

  integer :: i, info

  if (command_argument_count() >= 1) then
     call get_command_argument(1, solver_name, arg_len, status)
  else
     call rokko_serial_dense_ev_default_solver(solver_name)
  endif
  call rokko_split_solver_name(solver_name, library, routine)

  if (command_argument_count() == 2) then  
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     dim = 10
  endif

  print *,"library = ", library
  print *,"routine = ", routine
  print *,"dimension = ", dim

  call rokko_construct(solver, library)

  call rokko_construct(mat, dim, dim, rokko_matrix_col_major)
  call rokko_construct(Z, dim, dim, rokko_matrix_col_major)
  call rokko_construct(w, dim)

  ! generate frank matrix
  call rokko_frank_matrix_generate(mat)
  call rokko_print(mat)

  call rokko_construct(params)

  call rokko_set(params, "routine", "dsyevr")
  call rokko_set(params, "verbose", .true.)
  call rokko_set(params, "upper_index", 4)
  call rokko_set(params, "lower_index", 2)
  call rokko_set(params, "upper_value", 3.2d0)
  call rokko_set(params, "lower_value", 1.11d0)

!  call diagonalize(solver, mat, w, Z, params, params_out)
!  call diagonalize(solver, mat, w, Z, params)
!  call diagonalize(solver, mat, w, Z)
  call rokko_diagonalize(solver, mat, w, params, params_out)
!  call diagonalize(solver, mat, w, params)
!  call diagonalize(solver, mat, w)
  
  call rokko_get(params_out, "info", info)

  print *, "info=", info

  write(*,'(A)') "Computed Eigenvalues = "
  do i = 1, dim
     write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
  enddo

  call rokko_destruct(mat)
  call rokko_destruct(Z)
  call rokko_destruct(w)
  call rokko_destruct(solver)

end program frank_matrix
