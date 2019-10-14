!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

module mod_frank
  use iso_c_binding
  real(c_double), allocatable, dimension(:,:) :: eigen_array
  contains
    function func(i, j) bind(c)
      use iso_c_binding
      real(c_double) :: func
      integer(c_int), value, intent(in) :: i, j
      func = eigen_array(i+1,j+1)
    end function func
end module mod_frank

program frank_matrix_array_mpi
  use iso_c_binding
  use omp_lib
  use rokko
  use mod_frank
  use mpi
  implicit none
  integer :: dim
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_eigen_vector) :: w
  character(len=20) :: library, routine
  character(len=100) :: library_routine, tmp_str
  integer arg_len, status

  !---MPI variables---
  integer :: provided, ierr, myrank, nprocs

  !---loop variables---
  integer :: i, j
  integer :: m_global, n_global

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
  
  print *,"library = ", library
  print *,"routine = ", routine
  print *,"dimension = ", dim
  
  call rokko_parallel_dense_ev_construct(solver, library)
  call rokko_grid_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_parallel_dense_ev_default_mapping(solver, dim, grid, map)
  call rokko_distributed_matrix_construct(mat, map)
  call rokko_distributed_matrix_construct(Z, map)
  call rokko_eigen_vector_construct(w, dim)

  ! generate frank matrix
  m_global = rokko_distributed_matrix_get_m_global(mat)
  n_global = rokko_distributed_matrix_get_n_global(mat)

  allocate(eigen_array(dim,dim))
  do i = 1, m_global
     do j = 1, n_global
        eigen_array(j, i) = dim + 1 - max(i, j);
     end do
  end do
  call rokko_distributed_matrix_generate_function(mat, func)

  call rokko_distributed_matrix_print(mat)
  call rokko_parallel_dense_ev_diagonalize(solver, mat, w, Z)

  if (myrank.eq.0) then
     write(*,'(A)') "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_eigen_vector_get(w ,i)
     enddo
  endif

  call rokko_distributed_matrix_destruct(mat)
  call rokko_distributed_matrix_destruct(Z)
  call rokko_eigen_vector_destruct(w)
  call rokko_parallel_dense_ev_destruct(solver)
  call rokko_grid_destruct(grid)

  call MPI_finalize(ierr)

end program frank_matrix_array_mpi

