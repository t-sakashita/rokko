!/*****************************************************************************
!*
!* Rokko: Integrated Interface for libraries of eigenvalue decomposition
!*
!* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
!*
!* Distributed under the Boost Software License, Version 1.0. (See accompanying
!* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!*
!*****************************************************************************/

module mod_frank
  use iso_c_binding
  integer(c_int) :: dim_frank
  double precision, allocatable :: localized_array(:, :)
  contains
  double precision function func(i, j) bind(c)
    integer(c_int), value, intent(in) :: i, j
    print*, "i=", i
    func = localized_array(i+1,j+1)
  end function func
end module mod_frank

program frank_matrix
  use iso_c_binding
  use omp_lib
  use MPI
  use rokko
  use mod_frank
  implicit none

  integer::dim
  type(rokko_distributed_matrix) :: mat,Z !defined in rokko
  type(rokko_grid) :: grid !defined in rokko
  type(rokko_parallel_dense_solver) :: solver !defined in rokko

  type(rokko_localized_vector) :: w !localized_vector
  character(len=100) :: solver_name
  character(len=100) :: tmp_str
  integer args_cnt, arg_len, status

  !---MPI variables---
  integer :: provided,ierr,myrank,nprocs,comm,myrank_g,nprocs_g

  !---loop variables---
  integer :: i,j, count

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

  write(*,*) "solver_name=", solver_name
  write(*,*) "dim=",dim
  call rokko_parallel_dense_solver_construct(solver, solver_name)
  call rokko_grid_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)

  call rokko_distributed_matrix_construct(mat, dim, dim, grid, solver, rokko_matrix_col_major)
  call rokko_distributed_matrix_construct(Z, dim, dim, grid, solver, rokko_matrix_col_major)
  call rokko_localized_vector_construct(w, dim)

  ! generate frank matrix
  dim_frank = dim
  allocate(localized_array(dim,dim))
  do i = 1, dim
     do j = 1, dim
        localized_array(j, i) = dim + 1 - max(i, j);
     end do
  end do
  call rokko_distributed_matrix_generate_function_f(mat, func)

  call rokko_distributed_matrix_print(mat)
  call rokko_parallel_dense_solver_diagonalize_distributed_matrix(solver, mat, w, Z)

  if (myrank.eq.0) then
     write(*,'(A)') "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
     enddo
  endif

  call rokko_distributed_matrix_destruct(mat)
  call rokko_distributed_matrix_destruct(Z)
  call rokko_localized_vector_destruct(w)
  call rokko_parallel_dense_solver_destruct(solver)
  call rokko_grid_destruct(grid)

  call MPI_finalize(ierr)

end program frank_matrix
