!/*****************************************************************************
!*
!* Rokko: Integrated Interface for libraries of eigenvalue decomposition
!*
!* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!*                            Synge Todo <wistaria@comp-phys.org>,
!*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!*
!* Distributed under the Boost Software License, Version 1.0. (See accompanying
!* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!*
!*****************************************************************************/

program frank_matrix
  use MPI
  use rokko_sparse
  implicit none
  integer :: dim
  type(rokko_parallel_sparse_solver) :: solver
  type(rokko_distributed_crs_matrix) :: mat
  character(len=100) :: solver_name, tmp_str
  integer arg_len, status

  integer :: provided, ierr, myrank, nprocs
  integer :: i

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

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

  call rokko_parallel_sparse_solver_construct(solver, solver_name)

  call rokko_distributed_crs_matrix_construct(mat, dim, dim, solver)

  ! call rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(solver, mat)

  ! if (myrank.eq.0) then
  !    write(*,'(A)') "Computed Eigenvalues = "
  !    do i = 1, dim
  !       write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
  !    enddo
  ! endif

  call rokko_distributed_crs_matrix_destruct(mat)
  call rokko_parallel_sparse_solver_destruct(solver)

  call MPI_finalize(ierr)
end program frank_matrix
