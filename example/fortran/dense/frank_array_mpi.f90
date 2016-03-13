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
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_localized_vector) w
  character(len=100) :: solver_name, tmp_str
  integer args_cnt, arg_len, status
 
  double precision, allocatable, dimension(:,:) :: array

  integer :: provided,ierr, myrank, nprocs
  integer :: i, j

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

  call rokko_parallel_dense_ev_construct(solver, solver_name)
  call rokko_grid_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_parallel_dense_ev_default_mapping(solver, dim, grid, map)
  call rokko_distributed_matrix_construct(mat, map)
  call rokko_distributed_matrix_construct(Z, map)
  call rokko_localized_vector_construct(w, dim)

  ! generate frank matrix as a localized matrix
  allocate(array(dim, dim))
  do i=1, dim
     do j=1, dim
        array(i,j) = dble(dim + 1 - max(i, j))
     end do
  end do

!!$  do i=0, dim -1 
!!$    do j=0, dim -1
!!$      call rokko_distributed_matrix_set_global(mat,i,j,array(i+1,j+1))
!!$    enddo
!!$  enddo
  call rokko_distributed_matrix_generate_array(mat,array)
!  call generate_array_distributed_matrix(array, mat, dim, dim, dim)
  call rokko_distributed_matrix_print(mat)

  call rokko_parallel_dense_ev_diagonalize(solver, mat, w, Z)
!  call generate_distributed_matrix_array(Z, array, dim, dim, dim)

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_localized_vector_get(w, i)
     enddo
  endif

!  if (myrank.eq.0) then
!     print*, "array=", array
!  endif

  call rokko_distributed_matrix_destruct(mat)
  call rokko_distributed_matrix_destruct(Z)
  call rokko_localized_vector_destruct(w)
  call rokko_parallel_dense_ev_destruct(solver)
  call rokko_grid_destruct(grid)
  deallocate(array)

  call MPI_finalize(ierr)
end program frank_matrix
