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
module mod_frank
  integer :: n
  contains
  real*8 function func(i, j)
    integer, intent(in) :: i, j
!    print*, "i=", i
    func = n - max(dble(i), dble(j))
  end function func
end module mod_frank

real*8 function frank_matrix_element(i, j)
  integer, intent(in) :: i, j
  integer n
  common /frank_matrix_common/ n 
  frank_matrix_element = n - max(dble(i), dble(j))
end function frank_matrix_element

subroutine frank_matrix_set_dimension(n_in)
  integer n
  common /frank_matrix_common/ n
  n = n_in
end subroutine frank_matrix_set_dimension

program frank_matrix
  use MPI
  use rokko
  implicit none

  integer::dim
  type(distributed_matrix)::mat,Z !defined in rokko
  type(grid)::g !defined in rokko
  type(solver)::solver_ !defined in rokko
  
  real(8),allocatable::w(:),vec(:) !localized_vector
  character(len=100)::solver_name
  character(len=100)::tmp_str
  integer args_cnt, arg_len, status

  real*8 frank_matrix_element
  external frank_matrix_element

  !---MPI variables--
  integer::ierr,myrank,nprocs,comm,myrank_g,nprocs_g

  !---loop variables---
  integer::i,count

  call MPI_init(ierr) 
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  comm = MPI_COMM_WORLD

  args_cnt = command_argument_count ()
  if (args_cnt >= 2) then
     call get_command_argument (1, tmp_str, arg_len, status)
     solver_name = trim(tmp_str)
     call get_command_argument (2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     write(*,*) "error: command line"
     write(*,*) "usage: solver_name matrix_size"
     stop
  endif

  write(*,*) "solver_name=", solver_name  
  write(*,*) "dim=",dim
  call frank_matrix_set_dimension(dim)

  call set_solver(solver_, solver_name)
  write(*,*) "finished solver generation"
  call set_grid(g, MPI_COMM_WORLD)
  write(*,*) "finished grid generation"
  myrank_g = get_myrank_grid(g)
  nprocs_g = get_nprocs_grid(g)

  write(*,*) myrank_g,nprocs_g,myrank,nprocs

  call set_distributed_matrix(mat, dim, dim, g, solver_)
  call set_distributed_matrix(Z, dim, dim, g, solver_)
  allocate(w(dim));

  write(*,*) "finished matrix generation"
  
  call generate_distributed_matrix_function(mat, frank_matrix_element)

  call diagonalize(solver_, mat, w, Z)
  
  write(*,*) "finised matrix generation frank"

  if (myrank_g .eq. 0) then
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

