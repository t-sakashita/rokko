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

module various_mpi_comm_mod
  use MPI
  implicit none

contains

  integer function create_even_comm(comm_in)
    implicit none
    integer :: comm_in
    integer :: world_group, even_group
    integer :: nprocs, num_even
    integer, allocatable, dimension(:) :: even_members
    integer :: i, ierr

    call mpi_comm_size(comm_in, nprocs, ierr)

    num_even = (nprocs+1)/2
    allocate( even_members(num_even) )
    do i=0, num_even-1
       even_members(i+1) = 2 * i
    enddo
    call mpi_comm_group(MPI_COMM_WORLD, world_group, ierr)
    call mpi_group_incl(world_group, num_even, even_members, even_group, ierr)
    call mpi_comm_create(MPI_COMM_WORLD, even_group, create_even_comm, ierr)
    call mpi_group_free(world_group, ierr)
  end function create_even_comm

  integer function create_odd_comm(comm_in)
    implicit none
    integer :: comm_in
    integer :: world_group, odd_group
    integer :: nprocs, num_odd
    integer, allocatable, dimension(:) :: odd_members
    integer :: i, ierr

    call mpi_comm_size(comm_in, nprocs, ierr)

    num_odd = nprocs / 2
    allocate( odd_members(num_odd) )
    do i=0, num_odd-1
       odd_members(i+1) = 2 * i + 1
    enddo
    call mpi_comm_group(MPI_COMM_WORLD, world_group, ierr)
    call mpi_group_incl(world_group, num_odd, odd_members, odd_group, ierr)
    call mpi_comm_create(MPI_COMM_WORLD, odd_group, create_odd_comm, ierr)
    call mpi_group_free(world_group, ierr)
  end function create_odd_comm

  integer function create_even_odd_comm(comm_in) result(comm)
    implicit none
    integer :: comm_in
    integer :: myrank, ierr

    call mpi_comm_rank(comm_in, myrank, ierr)

    if (mod(myrank,2) == 0) then
       comm = create_even_comm(comm_in)
    else
       comm = create_odd_comm(comm_in)
    endif
  end function create_even_odd_comm

  integer function create_even_odd_comm_by_split(comm_in) result(comm)
    implicit none
    integer :: comm_in
    integer :: myrank, color
    integer :: ierr

    call mpi_comm_rank(comm_in, myrank, ierr)
    color = mod(myrank,2)
    call mpi_comm_split(comm_in, color, myrank, comm, ierr)
  end function create_even_odd_comm_by_split

  integer function find_square_root_like_divisor(n) result(i)
    integer, intent(in) :: n

    do i = int(sqrt(real(n))), 2, -1
       if (mod(n,i) == 0)  exit
    enddo
  end function find_square_root_like_divisor

  integer function create_cart_comm(comm_in) result(comm)
    implicit none
    integer :: comm_in
    integer :: nprocs
    integer :: dims(2)
    logical :: periods(2), reorder
    integer :: ierr

    call mpi_comm_size(comm_in, nprocs, ierr)
    dims(1) = find_square_root_like_divisor(nprocs)
    dims(2) = nprocs / dims(1)

    periods(1) = .false.
    periods(2) = .false.
    reorder = .false.
    call mpi_cart_create(comm_in, 2, dims, periods, reorder, &
         & comm, ierr)
  end function create_cart_comm

end module various_mpi_comm_mod
