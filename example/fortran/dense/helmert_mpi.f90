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

program helmert_matrix
  use rokko
  use mpi
  implicit none
  integer :: dim
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_eigen_vector) :: diag, w
  character(len=:), allocatable :: library, routine
  character(len=:), allocatable :: library_routine, tmp_str
  double precision, pointer, dimension(:) :: diag_ptr

  integer :: provided, ierr, myrank, nprocs
  integer :: i

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() >= 1) then
     call get_command_argument_deferred(1, library_routine)
  else
     library_routine = rokko_parallel_dense_ev_default_solver()
  endif
  call rokko_split_solver_name(library_routine, library, routine)

  if (command_argument_count() == 2) then
     call get_command_argument_deferred(2, tmp_str)
     read(tmp_str, *) dim
  else
     dim = 10
  endif

  if (myrank == 0) then
     print *,"library = ", library
     print *,"routine = ", routine
     print *,"dimension = ", dim
  endif

  call rokko_construct(solver, library)
  call rokko_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_default_mapping(solver, dim, grid, map)
  call rokko_construct(mat, map)
  call rokko_construct(Z, map)
  call rokko_construct(w, dim)

  ! generate matrix whose eigenvectors are helmert matrix
  call rokko_construct(diag, dim)
  call rokko_get_array_pointer(diag, diag_ptr)
  do i=1, dim
     diag_ptr(i) = i
  enddo
  call rokko_helmert_matrix_generate(mat, diag)
  call rokko_print(mat)

  call rokko_diagonalize(solver, mat, w, Z)

  if (myrank.eq.0) then
     write(*,'(A)') "Computed Eigenvalues = "
     call rokko_print(w)
     !     do i = 1, dim
     !        write(*,"(f30.20)") rokko_get_elem(w, i)
     !     enddo
  endif

  call rokko_destruct(mat)
  call rokko_destruct(Z)
  call rokko_destruct(w)
  call rokko_destruct(solver)
  call rokko_destruct(map)
  call rokko_destruct(grid)

  call MPI_finalize(ierr)
end program helmert_matrix
