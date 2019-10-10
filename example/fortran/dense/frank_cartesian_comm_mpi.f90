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
  integer, dimension(2) :: dims
  logical, dimension(2) :: periods 
  logical :: reorder
  
  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  dims(1) = int(sqrt(real(nprocs)));
  do while (.true.)
     if ( dims(1) == 1 ) exit
     if ( mod(nprocs, dims(1)) == 0 ) exit
     dims(1) = dims(1) - 1
  enddo
  dims(2) = nprocs / dims(1);
  periods(1) = .false.;  periods(2) = .false.
  reorder = .false.
  call mpi_cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, comm, ierr)
  if (myrank == 0) then
     write(*,'("Created ", i0, "x", i0, " size communicator with new cartesian topology")') dims(1), dims(2)
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

     if (myrank == 0) then
        print *,"library = ", library
        print *,"routine = ", routine
        print *,"dimension = ", dim
     endif
     
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
