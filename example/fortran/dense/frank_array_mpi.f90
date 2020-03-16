!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

program frank_matrix
  use rokko
  use mpi
  implicit none
  integer :: dim
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_eigen_vector) :: w
  character(len=:), allocatable :: library, routine
  character(len=:), allocatable :: library_routine, tmp_str
 
  double precision, allocatable, dimension(:,:) :: array
  integer :: provided,ierr, myrank, nprocs
  integer :: m_size, n_size, m_local, n_local
  integer :: local_i, local_j, global_i, global_j
  integer :: i

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() >= 1) then
     call get_command_argument_deferred(1, library_routine)
  else
     call rokko_parallel_dense_ev_default_solver(library_routine)
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

  call rokko_diagonalize(solver, mat, w, Z)

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_get_elem_f(w, i)
     enddo
  endif

  deallocate(array)
  call rokko_destruct(mat)
  call rokko_destruct(Z)
  call rokko_destruct(w)
  call rokko_destruct(solver)
  call rokko_destruct(grid)

  call MPI_finalize(ierr)
end program frank_matrix
