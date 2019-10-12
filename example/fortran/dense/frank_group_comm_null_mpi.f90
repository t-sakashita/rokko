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

program frank_matrix
  use rokko
  implicit none
  include 'mpif.h'
  integer :: dim
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_eigen_vector) :: w
  character(len=20) :: library, routine
  character(len=100) :: library_routine, tmp_str
  integer arg_len, status

  integer :: provided, ierr, myrank, nprocs
  integer :: i
  integer :: comm
  integer :: world_group, even_group
  integer :: num_even
  integer, allocatable,dimension(:) :: even_members
    
  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  num_even = (nprocs+1)/2
  allocate( even_members(num_even) )
  do i=0, num_even-1
     even_members(i+1) = 2 * i     
  enddo
  call mpi_comm_group(MPI_COMM_WORLD, world_group, ierr)
  call mpi_group_incl(world_group, num_even, even_members, even_group, ierr)
  call mpi_comm_create(MPI_COMM_WORLD, even_group, comm, ierr)
  call mpi_group_free(world_group, ierr)
  if (comm == MPI_COMM_NULL) then
     print *,"orig_rank=", myrank, " is COMM_NULL"
  else
     print *,"orig_rank=", myrank, " is NOT COMM_NULL"
  endif

  if (comm /= MPI_COMM_NULL) then
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
     
     call rokko_construct(solver, library)
     call rokko_construct(grid, comm, rokko_grid_row_major)
     call rokko_default_mapping(solver, dim, grid, map)
     call rokko_construct(mat, map)
     call rokko_construct(Z, map)
     call rokko_construct(w, dim)
     
     ! generate frank matrix
     call rokko_frank_matrix_generate(mat)
     call rokko_print(mat)
     
     call rokko_diagonalize(solver, mat, w, Z)
     
     if (myrank.eq.0) then
        write(*,'(A)') "Computed Eigenvalues = "
        do i = 1, dim
           write(*,"(f30.20)") rokko_eigen_vector_get(w, i)
        enddo
     endif
     
     call rokko_destruct(mat)
     call rokko_destruct(Z)
     call rokko_destruct(w)
     call rokko_destruct(solver)
     call rokko_destruct(map)
     call rokko_destruct(grid)
     call mpi_comm_free(comm, ierr)
  endif
  call MPI_finalize(ierr)
end program frank_matrix
   
