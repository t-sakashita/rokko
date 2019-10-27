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

module heisenberg
  use mpi
  use rokko
  implicit none
  integer :: comm
  integer :: nprocs, myrank
  integer, private :: dim, local_offset, num_local_rows
  double precision, allocatable, private :: buffer(:)
  integer, private :: L, lattice_size, p
  integer, allocatable, private :: lattice(:,:)

contains

  subroutine initialize (mat, L_in, lattice_size_in, lattice_in)
    type(rokko_distributed_mfree), intent(inout) :: mat
    integer, intent(in) :: L_in, lattice_size_in
    integer :: lattice_in(2,lattice_size_in)
    integer :: ierr

    L = L_in
    lattice_size = lattice_size_in

    comm = mpi_comm_world
    call mpi_comm_rank(comm, myrank, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    p = find_power_of_two(nprocs)
    if (nprocs /= (2 ** p)) then
       if ( myrank == 0 ) then
          write(0,*) "Error: This program can be run only for powers of 2"
       endif
       call mpi_abort(comm, 1, ierr)
    endif

    dim = 2 ** L
    num_local_rows = 2 ** (L-p)
    local_offset = num_local_rows * myrank
    !print*, "myrank=", myrank, "num_local_rows=", num_local_rows

    allocate( lattice(2,lattice_size) )
    lattice = lattice_in
    allocate( buffer(0:num_local_rows-1) )
    call rokko_distributed_mfree_construct(mat, multiply, dim, num_local_rows)
  end subroutine initialize

  integer function find_power_of_two(n_in)
    integer, intent(in) :: n_in
    integer :: n

    find_power_of_two = 0
    n = n_in
    do while (mod(n,2) == 0)
       n = n / 2
       find_power_of_two = find_power_of_two + 1
    enddo
  end function find_power_of_two

  integer function get_num_local_rows ()
    get_num_local_rows = num_local_rows
  end function get_num_local_rows

  integer function get_dim ()
    get_dim = dim
  end function get_dim

  ! subroutine passed to c function.
  ! it must be interoperable!
  subroutine multiply (n, x, y)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(0:n-1)
    double precision, intent(out) :: y(0:n-1)
    integer :: ierr
    integer :: status(mpi_status_size)
    integer :: lattice_i, i, j, m1, m2, m3, m
    integer :: k

    if (num_local_rows == 0) then
       return
    endif

    do lattice_i = 1, lattice_size
       i = lattice(1,lattice_i)
       j = lattice(2,lattice_i)
       if (i < (L-p)) then
          if (j < (L-p)) then
             m1 = 2 ** i
             m2 = 2 ** j
             m3 = m1 + m2
#ifdef _OPENMP
!$omp do
#endif
             do k=0, n-1
                if (and(k,m3) == m1) then  ! when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
                   y(k) = y(k) + 0.5 * x(xor(k,m3)) - 0.25 * x(k)
                else if (and(k,m3) == m2) then
                   y(k) = y(k) + 0.5 * x(xor(k,m3)) - 0.25 * x(k)
                else
                   y(k) = y(k) + 0.25 * x(k)
                endif
             enddo
          else
             m = 2 ** (j-(L-p))
             call MPI_Sendrecv(x, N, MPI_DOUBLE, xor(myrank,m), 0, &
                  &            buffer, N, MPI_DOUBLE, xor(myrank,m), 0, &
                  &            comm, status, ierr)
             m1 = 2 ** i
             if (and(myrank,m) == m) then
#ifdef _OPENMP
!$omp do
#endif
                do k=0, n-1
                   if (and(k,m1) == m1) then
                      y(k) = y(k) + 0.25 * x(k)
                   else
                      y(k) = y(k) + 0.5 * buffer(xor(k,m1)) - 0.25 * x(k)
                   endif
                enddo
             else
#ifdef _OPENMP
!$omp do
#endif
                do k=0, n-1
                   if (and(k,m1) == m1) then
                      y(k) = y(k) + 0.5 * buffer(xor(k,m1)) - 0.25 * x(k)
                   else
                      y(k) = y(k) + 0.25 * x(k)
                   endif
                enddo
             endif
          endif
       else
          if (j < (L-p)) then
             m = 2 ** (i-(L-p))
             call MPI_Sendrecv(x, N, MPI_DOUBLE, xor(myrank,m), 0, &
                  &            buffer, N, MPI_DOUBLE, xor(myrank,m), 0, &
                  &            comm, status, ierr)
             m1 = 2 ** j
             if (and(myrank,m) == m) then
#ifdef _OPENMP
!$omp do
#endif
                do k=0, n-1
                   if (and(k,m1) == m1) then
                      y(k) = y(k) + 0.25 * x(k)
                   else
                      y(k) = y(k) + 0.5 * buffer(xor(k,m1)) - 0.25 * x(k)
                   endif
                enddo
             else
#ifdef _OPENMP
!$omp do
#endif
                do k=0, n-1
                   if (and(k,m1) == m1) then
                      y(k) = y(k) + 0.5 * buffer(xor(k,m1)) - 0.25 * x(k)
                   else
                      y(k) = y(k) + 0.25 * x(k)
                   endif
                enddo
             endif
          else
             m = (2 ** (i-(L-p))) + (2 ** (j-(L-p)))
             if ((and(myrank,m) /= m) .and. (and(myrank,m) /= 0)) then
                call MPI_Sendrecv(x, N, MPI_DOUBLE, xor(myrank,m), 0, &
                     &            buffer, N, MPI_DOUBLE, xor(myrank,m), 0, &
                     &            comm, status, ierr)
#ifdef _OPENMP
!$omp do
#endif
                do k=0, n-1
                   y(k) = y(k) + 0.5 * buffer(k) - 0.25 * x(k)
                enddo
             else
#ifdef _OPENMP
!$omp do
#endif
                do k=0, n-1
                   y(k) = y(k) + 0.25 * x(k)
                enddo
             endif
          endif
       endif
    enddo
  end subroutine multiply

end module heisenberg

program main
  use mpi
  use rokko
  use heisenberg
  implicit none
  integer :: provided, ierr

  type(rokko_parameters) :: params, params_out
  double precision :: eig_val
  double precision, allocatable, dimension(:) :: eig_vec

  type(rokko_parallel_sparse_ev) :: solver
  character(len=100) :: library_routine, tmp_str
  character(len=50) :: library, routine
  integer arg_len, status
  type(rokko_distributed_mfree) :: mat
  integer :: dim, i
  integer :: num_evals, block_size, max_iters
  integer :: num_local_rows, num_conv
  integer :: L, lattice_size
  integer, allocatable :: lattice(:,:)
  double precision :: tol

  call mpi_init_thread(mpi_thread_multiple, provided, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)

  if (command_argument_count() >= 1) then
     call get_command_argument(1, library_routine, arg_len, status)
  else
     call rokko_parallel_sparse_ev_default_solver(library_routine)
  endif
  call rokko_split_solver_name(library_routine, library, routine)

  if (command_argument_count() == 2) then
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) L
  else
     L = 10  ! default
  endif

  lattice_size = L
  allocate( lattice(2,lattice_size) )
  do i=0, L-1
     lattice(1,i+1) = i
     lattice(2,i+1) = mod(i+1, L)
  end do

  call rokko_construct(solver, library)
  call initialize(mat, L, L, lattice)
  dim = get_dim()
  num_local_rows = get_num_local_rows()
  allocate( eig_vec(num_local_rows) )

  if (myrank == 0) then
     write(*,*) "solver name = ", trim(library)
     write(*,*) "matrix dimension = ", dim
  endif

  call rokko_construct(params)
  call rokko_set(params, "routine", routine)
  call rokko_set(params, "verbose", .true.)
  call rokko_set(params, "num_evals", 1)
  call rokko_set(params, "block_size", 5)
  call rokko_set(params, "max_iters", 500)
  call rokko_set(params, "conv_tol", 1.0d-3)
  call rokko_set(params, "verbose", .true.)
  call rokko_diagonalize(solver, mat, params, params_out)

  num_conv = rokko_parallel_sparse_ev_num_conv(solver)
  if (num_conv == 0) then
     call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  endif

  if (myrank == 0) then
     eig_val = rokko_parallel_sparse_ev_eigenvalue(solver, 0)
     call rokko_parallel_sparse_ev_eigenvector(solver, 0, eig_vec)
     print *, "eigval=", eig_val
     print *, "Computed Eigenvector = "
     print *, eig_vec
  endif
  call rokko_destruct(mat)
  call rokko_destruct(solver)

  call mpi_finalize(ierr)
end program main
