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

program generate_matrix012_mpi
  use rokko
  use mpi
  implicit none
  integer :: dim
  type(rokko_distributed_matrix) :: mat
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map

  integer :: provided, ierr, myrank, nprocs

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  dim = 10

  if (myrank == 0) then
     print *,"dimension = ", dim
  endif

  call rokko_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_construct(map, dim, 1, grid)
  call rokko_construct(mat, map)

  ! generate matrix012
  call rokko_matrix012_generate(mat)
  call rokko_print(mat)

  call rokko_destruct(mat)
  call rokko_destruct(map)
  call rokko_destruct(grid)

  call MPI_finalize(ierr)
end program generate_matrix012_mpi
