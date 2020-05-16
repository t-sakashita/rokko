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

program main
  use mpi
  implicit none
  integer, parameter :: n = 8, nb = 1
  integer :: provided, ierr, myrank, nprocs
  integer :: icontxt, npcol, nprow, desc(9)
  integer :: i, j, info
  real(4), allocatable :: a(:, :), z(:,:), w(:)
  integer :: lwork, liwork
  real(4), allocatable :: work(:)
  integer, allocatable :: iwork(:)

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  nprow = int(sqrt(1. * nprocs + 0.5))
  npcol = int(1. * nprocs / nprow)
  if (myrank == 0) then
     print *, "n =", n
     print *, "nprocs =", nprocs
     print *, "nprow =", nprow
     print *, "npcol =", npcol
     if ((nprocs /= nprow * npcol) .or. (mod(n, nprow) /= 0) .or. &
          (mod(n, npcol) /= 0)) then
        print *, "incompatible matrix size and number of processes"
        call MPI_Abort(MPI_COMM_WORLD, 127, ierr)
     end if
  end if
  call BLACS_get(0, 0, icontxt)
  call BLACS_gridinit(icontxt, 'R', nprow, npcol)

  call descinit(desc, n, n, nb, nb, 0, 0, icontxt, n/nprow, info)
  allocate(a(n/nprow, n/npcol), z(n/nprow, n/npcol), w(n))
  do j = 1, n
     do i = 1, n
        call pselset(a, i, j, desc, real(min(i, j)))
     end do
  end do

  lwork = -1
  liwork = 1
  allocate(work(1), iwork(1))
  call pssyevd('V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, work, lwork, iwork, liwork, info)
  lwork = int(work(1))
  liwork = iwork(1)
  deallocate(work, iwork)
  allocate(work(lwork), iwork(liwork))

  call pssyevd('V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, work, lwork, iwork, liwork, info)
  if (myrank == 0) then
     print *, "eigenvalues:", w(:)
  end if

  deallocate(a, z, w, work, iwork)
  call BLACS_gridexit(icontxt)
  call MPI_finalize(ierr)
end program main
