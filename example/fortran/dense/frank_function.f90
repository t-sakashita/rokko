!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

module frank_mod
  implicit none
  public frank_matrix_element
  integer, private :: dim

contains

  function frank_matrix_element(i, j)
    double precision :: frank_matrix_element
    integer, intent(in) :: i, j
    frank_matrix_element = dble(dim - max(i, j))
  end function frank_matrix_element

  subroutine frank_matrix_set_dimension(dim_in)
    integer, intent(in) :: dim_in
    dim = dim_in
  end subroutine frank_matrix_set_dimension

end module frank_mod

program frank_matrix
  use rokko
  use frank_mod
  implicit none
  integer :: dim
  type(rokko_serial_dense_ev) :: solver
  type(rokko_eigen_matrix) :: mat, Z
  type(rokko_eigen_vector) :: w
  type(rokko_parameters) :: params, params_out
  character(len=:), allocatable :: library, routine
  character(len=:), allocatable :: solver_name, tmp_str

  integer :: i, info

  if (command_argument_count() >= 1) then
     call get_command_argument_deferred(1, solver_name)
  else
     solver_name = rokko_serial_dense_ev_default_solver()
  endif
  call rokko_split_solver_name(solver_name, library, routine)

  if (command_argument_count() == 2) then
     call get_command_argument_deferred(2, tmp_str)
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
  call frank_matrix_set_dimension(dim)
  call rokko_generate0(mat, frank_matrix_element)
  call rokko_print(mat)

  call rokko_construct(params)
  call rokko_set(params, "routine", "dsyev")

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
     write(*,"(f30.20)") rokko_get_elem(w ,i)
  enddo

  call rokko_destruct(mat)
  call rokko_destruct(Z)
  call rokko_destruct(w)
  call rokko_destruct(solver)

end program frank_matrix
