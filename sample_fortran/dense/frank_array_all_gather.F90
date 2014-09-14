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
  use rokko
  implicit none
  integer :: dim
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_grid) :: grid
  type(rokko_parallel_dense_solver) :: solver
  type(rokko_localized_vector) :: w
  character(len=100) :: solver_name, tmp_str
  integer arg_len, status
  real(8), allocatable, dimension(:,:) :: array, array_tmp

  integer :: provided,ierr, myrank, nprocs
  integer :: i, j, proc

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count().eq.2) then
     call get_command_argument(1, tmp_str, arg_len, status)
     solver_name = trim(tmp_str)
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     write(*,'(A)') "Error: frank_array solver_name dimension"
     stop
  endif

  write(*,*) "solver name = ", trim(solver_name)
  write(*,*) "matrix dimension = ", dim

  call rokko_parallel_dense_solver_construct(solver, solver_name)
  call rokko_grid_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)

  call rokko_distributed_matrix_construct(mat, dim, dim, grid, solver, rokko_matrix_col_major)
  call rokko_distributed_matrix_construct(Z, dim, dim, grid, solver, rokko_matrix_col_major)
  call rokko_localized_vector_construct(w, dim)

  ! generate frank matrix as a localized matrix
  allocate(array(dim, dim))
  do i=1, dim
     do j=1, dim
        array(i,j) = dble(dim + 1 - max(i, j))
     end do
  end do
  allocate(array_tmp(dim, dim))

  call rokko_distributed_matrix_generate_array(mat,array)
  call rokko_distributed_matrix_print(mat)
  call rokko_parallel_dense_solver_diagonalize_distributed_matrix(solver, mat, w, Z)

  array = 0.0
  call rokko_all_gather(Z, array_tmp)
  array = matmul(transpose(array_tmp), array_tmp)
  call mpi_barrier(mpi_comm_world, ierr)                                                                                                                                       

  !  print*, "array=", array
  do proc = 0, nprocs-1
     if (proc == myrank) then
        print*, "fmyrank=", myrank
     do i = 1, dim
        write(*,'(10f8.4)') (array(i, j), j=1, dim)
     end do
     endif
   call mpi_barrier(mpi_comm_world, ierr)
!   call sleep(0.1)
  end do

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
     enddo
  endif

  call mpi_barrier(mpi_comm_world, ierr)                                                                                                                                       

  call rokko_distributed_matrix_destruct(mat)
  call rokko_distributed_matrix_destruct(Z)
  call rokko_localized_vector_destruct(w)
  call rokko_parallel_dense_solver_destruct(solver)
  call rokko_grid_destruct(grid)
  deallocate(array,array_tmp)

  call MPI_finalize(ierr)
end program frank_matrix
