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
  use rokko_f
  implicit none
  integer :: dim
  type(distributed_matrix) :: mat, Z
  type(grid) :: g
  type(solver) :: solver_
  real(8), allocatable :: w(:), vec(:)
  character(len=100) :: solver_name, tmp_str
  integer args_cnt, arg_len, status
  real(8), allocatable, dimension(:,:) :: array

  integer :: ierr, myrank, nprocs
  integer :: i, j

  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count().eq.2) then
     call get_command_argument(1, tmp_str, arg_len, status)
     solver_name = trim(tmp_str)
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     write(*,'(A)') "Error: frank_array solver_name dimension"
     stop
  endif

  write(*,*) "solver name = ", trim(solver_name)
  write(*,*) "matrix dimension = ", dim

  call set_solver(solver_, solver_name)
  call set_grid(g, MPI_COMM_WORLD)

  call set_distributed_matrix(mat, dim, dim, g, solver_)
  call set_distributed_matrix(Z, dim, dim, g, solver_)
  allocate(w(dim));

  ! generate frank matrix as a localized matrix
  allocate(array(dim, dim))
  do i=1, dim
     do j=1, dim
        array(i,j) = dble(dim + 1 - max(i, j))
     end do
  end do

  call generate_array_distributed_matrix(array, mat, dim, dim, dim)
  call print_distributed_matrix(mat)

  call diagonalize(solver_, mat, w, Z)
!  call generate_distributed_matrix_array(Z, array, dim, dim, dim)

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") w(i)
     enddo
  endif

!  if (myrank.eq.0) then
!     print*, "array=", array
!  endif

  call del_distributed_matrix(mat)
  call del_distributed_matrix(Z)
  call del_solver(solver_)
  deallocate(w)
  deallocate(array)

  call MPI_finalize(ierr)
end program frank_matrix
