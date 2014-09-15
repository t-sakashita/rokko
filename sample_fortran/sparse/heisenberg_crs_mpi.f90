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
  use rokko_sparse
  implicit none
  integer :: provided, ierr, myrank, nprocs

  integer :: L, k
  integer :: dim
  integer, allocatable, dimension(:) :: lattice_first, lattice_second
  real(8) :: diag
  integer :: i, j, m1, m2, m3, count
  integer :: row, start_row, end_row
  integer(c_int), allocatable, target, dimension(:) :: cols
  real(c_double), allocatable, target, dimension(:) :: values
!  integer(4), allocatable, dimension(:) :: cols
!  real(8), allocatable, dimension(:) :: values

  type(rokko_parallel_sparse_solver) :: solver
  type(rokko_distributed_crs_matrix) :: mat
  character(len=100) :: solver_name, tmp_str
  integer :: arg_len, status
  integer :: num_evals, block_size, max_iters
  real(8) :: tol

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  solver_name = "anasazi"

  L = 3
  dim = ishft(1,L)
  allocate( lattice_first(L) )
  allocate( lattice_second(L) )
  do k = 1, L
     lattice_first(k) = k - 1
     lattice_second(k) = mod(k, L)
  end do

  write(*,*) "solver name = ", trim(solver_name)
  write(*,*) "matrix dimension = ", dim

  call rokko_parallel_sparse_solver_construct(solver, solver_name)

  call rokko_distributed_crs_matrix_construct(mat, dim, dim, solver)

  start_row = rokko_distributed_crs_matrix_start_row(mat);
  end_row = rokko_distributed_crs_matrix_end_row(mat);
  print*, "row_start=", start_row
  print*, "row_end=", end_row

  allocate( cols(dim) )
  allocate( values(dim) )

  do row = start_row, end_row
     count = 0
     diag = 0;
     do k = 1, L
        i = lattice_first(k)
        j = lattice_second(k)
        m1 = ishft(1, i)
        m2 = ishft(1, j)
        m3 = m1 + m2
        if ( iand(row, m3) == m1 .or. iand(row, m3) == m2 ) then
           cols(count) = ieor(row, m3)
           values(count) = 0.5
           count = count + 1
           diag = diag - 0.25
        else
           diag = diag + 0.25
        end if
     end do
     cols(count) = row
     values(count) = diag
     count = count + 1
     print*, "cols=", cols
     print*, "values=", values
     call rokko_distributed_crs_matrix_insert(mat, row, count, cols, values);
  end do
  call rokko_distributed_crs_matrix_complete(mat);

  num_evals = 10;
  block_size = 5;
  max_iters = 500;
  tol = 1.0e-8;
  call rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(solver, mat, num_evals, block_size, max_iters, tol)

  ! if (myrank.eq.0) then
  !    write(*,'(A)') "Computed Eigenvalues = "
  !    do i = 1, dim
  !       write(*,"(f30.20)") rokko_localized_vector_get(w ,i)
  !    enddo
  ! endif

  call rokko_distributed_crs_matrix_destruct(mat)
  call rokko_parallel_sparse_solver_destruct(solver)

  call MPI_finalize(ierr)
end program frank_matrix
