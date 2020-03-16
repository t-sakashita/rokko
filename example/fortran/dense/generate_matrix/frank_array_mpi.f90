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

program generate_frank_array_mpi
  use rokko
  use mpi
  implicit none
  integer :: dim
  type(rokko_distributed_matrix) :: mat
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  double precision, allocatable :: array(:,:)
  integer :: provided, ierr, myrank, nprocs
  integer :: m_size, n_size, m_local, n_local
  integer :: local_i, local_j, global_i, global_j

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  dim = 10
  if (myrank == 0) then
     print *,"dimension = ", dim
  endif

  call rokko_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_construct(map, dim, 1, grid)
  m_size = rokko_get_m_size(map)
  n_size = rokko_get_n_size(map)
  allocate(array(m_size, n_size))

  ! generate frank matrix
  m_local = rokko_get_m_local(map)
  n_local = rokko_get_n_local(map)
  do local_j = 0, n_local-1
     global_j = rokko_translate_l2g_col0(map, local_j)
     do local_i = 0, m_local-1
        global_i = rokko_translate_l2g_row0(map, local_i)
        array(local_i+1, local_j+1) = dim - max(global_i, global_j)
     enddo
  enddo

  call rokko_construct(mat, map, array)
  call rokko_print(mat)

  deallocate(array)
  call rokko_destruct(mat)
  call rokko_destruct(map)
  call rokko_destruct(grid)

  call MPI_finalize(ierr)
end program generate_frank_array_mpi
