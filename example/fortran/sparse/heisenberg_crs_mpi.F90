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

program heisenberg_crs_mpi
  use rokko
  use mpi
  implicit none
  integer :: provided, ierr, myrank, nprocs
  integer :: L, k
  integer :: dim
  integer, allocatable, dimension(:) :: lattice_first, lattice_second
  double precision :: diag
  integer :: i, j, m1, m2, m3, count
  integer :: row, start_row, end_row
  integer, allocatable, dimension(:) :: cols
  double precision, allocatable, dimension(:) :: values

  double precision :: eig_val
  double precision, allocatable, dimension(:) :: eig_vec

  type(rokko_parallel_sparse_ev) :: solver
  type(rokko_mapping_1d) :: map
  type(rokko_distributed_crs_matrix) :: mat
  character(len=:), allocatable :: library, routine
  character(len=:), allocatable :: library_routine, tmp_str
  type(rokko_parameters) :: params, params_out
  integer :: num_local_rows, num_conv

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() >= 1) then
     call get_command_argument_deferred(1, library_routine)
  else
     call rokko_parallel_sparse_ev_default_solver(library_routine)
  endif
  call rokko_split_solver_name(library_routine, library, routine)
  
  if (command_argument_count() == 2) then  
     call get_command_argument_deferred(2, tmp_str)
     read(tmp_str, *) L
  else
     L = 8
  endif
  
  dim = ishft(1, L)
  allocate( lattice_first(L) )
  allocate( lattice_second(L) )
  do k = 1, L
     lattice_first(k) = k - 1
     lattice_second(k) = mod(k, L)
  end do

  if (myrank == 0) then
     write(*,*) "solver name = ", trim(library)
     write(*,*) "matrix dimension = ", dim
  endif

  call rokko_construct(solver, library)

  call rokko_default_mapping(solver, dim, mpi_comm_world, map)
  call rokko_construct(mat, map, L+1)

  start_row = rokko_start_row0(map)  ! C style index
  end_row = rokko_end_row(map) - 1  ! Fortran style index

  allocate( cols(dim) )
  allocate( values(dim) )

  do row = start_row, end_row
     count = 0
     diag = 0
     do k = 1, L
        i = lattice_first(k)
        j = lattice_second(k)
        m1 = ishft(1, i)
        m2 = ishft(1, j)
        m3 = m1 + m2
        if ( (iand(row, m3) == m1) .or. (iand(row, m3) == m2) ) then
           count = count + 1
           cols(count) = ieor(row, m3)
           values(count) = 0.5
           diag = diag - 0.25
        else
           diag = diag + 0.25
        end if
     end do
     count = count + 1
     cols(count) = row
     values(count) = diag
     call rokko_insert0(mat, row, count, cols, values)
  end do
  call rokko_complete(mat)
!  call rokko_print(mat)

  call rokko_construct(params)
  call rokko_set(params, "routine", routine)
  call rokko_diagonalize(solver, mat, params)
  
  num_conv = rokko_num_conv(solver)
  if ((num_conv >= 1) .and. (myrank == 0)) then
     eig_val = rokko_eigenvalue(solver, 0)
     print *, "eigval=", eig_val
  endif

  eig_val = rokko_eigenvalue(solver, 0)
  num_local_rows = rokko_num_local_rows(mat)
  print*, "num_local_rows=", num_local_rows
  allocate( eig_vec(num_local_rows) )
  call rokko_eigenvector(solver, 0, eig_vec)

  if (myrank.eq.0) then
     print*, "number of converged eigenpairs=", num_conv
     print*, "Computed Eigenvalues = ", eig_val
     print*, "Computed Eigenvector = "
     print '(8f10.4)', eig_vec
  endif

  call rokko_destruct(params)
  call rokko_destruct(mat)
  call rokko_destruct(solver)

  call MPI_finalize(ierr)
end program heisenberg_crs_mpi

