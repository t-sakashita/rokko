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

real*8 function frank_matrix_element(i, j)
  integer, intent(in) :: i, j
  integer dim
  common /frank_matrix_common/ dim
  frank_matrix_element = dble(dim + 1 - max(i, j))
end function frank_matrix_element

subroutine frank_matrix_set_dimension(dim_in)
  integer, intent(in) :: dim_in
  integer dim
  common /frank_matrix_common/ dim
  dim = dim_in
end subroutine frank_matrix_set_dimension

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

  integer :: ierr, myrank, nprocs
  integer :: i

  real*8 frank_matrix_element
  external frank_matrix_element

  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count().eq.2) then
     call get_command_argument(1, tmp_str, arg_len, status)
     solver_name = trim(tmp_str)
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     write(*,'(A)') "Error: frank_function solver_name dimension"
     stop
  endif

  write(*,*) "solver name = ", trim(solver_name)
  write(*,*) "matrix dimension = ", dim

  call set_solver(solver_, solver_name)
  call set_grid(g, MPI_COMM_WORLD)

  call set_distributed_matrix(mat, dim, dim, g, solver_)
  call set_distributed_matrix(Z, dim, dim, g, solver_)
  allocate(w(dim));

  ! generate frank matrix from frank_matrix_element function
  call frank_matrix_set_dimension(dim)
  call generate_distributed_matrix_function(mat, frank_matrix_element)
  call print_distributed_matrix(mat)

  call diagonalize(solver_, mat, w, Z)

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") w(i)
     enddo
  endif

  call del_distributed_matrix(mat)
  call del_distributed_matrix(Z)
  call del_solver(solver_)
  deallocate(w)

  call MPI_finalize(ierr)
end program frank_matrix
