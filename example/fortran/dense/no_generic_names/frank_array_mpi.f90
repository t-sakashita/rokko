!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

program frank_matrix
  use rokko
  implicit none
  include 'mpif.h'
  integer :: dim
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_eigen_vector) w
  character(len=20) :: library, routine
  character(len=100) :: library_routine, tmp_str
  integer arg_len, status
 
  double precision, allocatable, dimension(:,:) :: array

  integer :: provided,ierr, myrank, nprocs
  integer :: i, j

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() >= 1) then
     call get_command_argument(1, library_routine, arg_len, status)
  else
     call rokko_parallel_dense_ev_default_solver(library_routine)
  endif
  call rokko_split_solver_name(library_routine, library, routine)

  if (command_argument_count() == 2) then  
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     dim = 10
  endif
  
  write(*,*) "solver name = ", trim(library)
  write(*,*) "matrix dimension = ", dim
  
  call rokko_parallel_dense_ev_construct(solver, library)
  call rokko_grid_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_parallel_dense_ev_default_mapping(solver, dim, grid, map)
  call rokko_distributed_matrix_construct(mat, map)
  call rokko_distributed_matrix_construct(Z, map)
  call rokko_eigen_vector_construct(w, dim)

  ! generate frank matrix as a eigen matrix
  allocate(array(dim, dim))
  do i=1, dim
     do j=1, dim
        array(i,j) = dble(dim + 1 - max(i, j))
     end do
  end do

  call rokko_distributed_matrix_generate_from_array(mat,array)
  call rokko_distributed_matrix_print(mat)

  call rokko_parallel_dense_ev_diagonalize(solver, mat, w, Z)

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_eigen_vector_get(w, i)
     enddo
  endif

!  if (myrank.eq.0) then
!     print*, "array=", array
!  endif

  call rokko_distributed_matrix_destruct(mat)
  call rokko_distributed_matrix_destruct(Z)
  call rokko_eigen_vector_destruct(w)
  call rokko_parallel_dense_ev_destruct(solver)
  call rokko_grid_destruct(grid)
  deallocate(array)

  call MPI_finalize(ierr)
end program frank_matrix
