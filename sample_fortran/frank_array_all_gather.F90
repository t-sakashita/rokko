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
  implicit none
  integer :: dim
  type(distributed_matrix) :: mat, Z
  type(grid) :: g
  type(solver) :: solver_
  real(8), allocatable :: w(:), vec(:)
  character(len=100) :: solver_name, tmp_str
  integer args_cnt, arg_len, status
  real(8), allocatable, dimension(:,:) :: array, array_tmp

  integer :: ierr, myrank, nprocs
  integer :: i, j, proc

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
  allocate(array_tmp(dim, dim))

  call generate_array_distributed_matrix(array, mat, dim, dim, dim)
!  call print_distributed_matrix(mat)
  call diagonalize(solver_, mat, w, Z)

  array = 0.0
  call all_gather(Z, array_tmp)
  array = matmul(transpose(array_tmp), array_tmp)
  call mpi_barrier(mpi_comm_world, ierr)                                                                                                                                       

  !  print*, "array=", array
  do proc = 0, nprocs-1
     if (proc == myrank) then
      print*, "fmyrank=", myrank
     do i = 1, dim
        write(*,'(10f8.4)') (array(i, j), j=1, dim)
     end do
     endif
   call mpi_barrier(mpi_comm_world, ierr)
!   call sleep(0.1)
  end do

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") w(i)
     enddo
  endif

!  if (myrank.eq.0) then
!     print*, "array=", array
!  endif
  write(*,*) "myrank:", myrank, " fortran_file:",__FILE__, " line:",&
             __LINE__
  call mpi_barrier(mpi_comm_world, ierr)                                                                                                                                       

  call del_distributed_matrix(mat)
  call del_distributed_matrix(Z)
  call del_solver(solver_)
  deallocate(w)
  deallocate(array)

  call MPI_finalize(ierr)
end program frank_matrix
