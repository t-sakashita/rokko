!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program frank_matrix
  use rokko
  use various_mpi_comm_mod
  implicit none
  integer :: dim
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_eigen_vector) :: w
  character(len=:), allocatable :: library, routine
  character(len=:), allocatable :: library_routine, tmp_str
  integer :: comm
  integer :: provided, ierr, myrank, nprocs
  integer :: i

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  comm = create_even_comm(mpi_comm_world)
  if (comm == MPI_COMM_NULL) then
     print *,"orig_rank=", myrank, " is COMM_NULL"
  else
     print *,"orig_rank=", myrank, " is NOT COMM_NULL"
  endif

  if (comm /= MPI_COMM_NULL) then
     if (command_argument_count() >= 1) then
        call get_command_argument_deferred(1, library_routine)
     else
        library_routine = rokko_parallel_dense_ev_default_solver()
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
           write(*,"(f30.20)") rokko_get_elem(w, i)
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
   
