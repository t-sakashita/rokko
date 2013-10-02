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
  use iso_c_binding
  contains
  real(c_double) function func(i, j) bind(c)
    integer(c_int) :: n
    common /mydata/n
    integer(c_int), value, intent(in) :: i, j
!    print*, "i=", i
    func = n - max(dble(i), dble(j))
  end function func
end module mod_frank

program frank_matrix
  use iso_c_binding
  use omp_lib
  use MPI
  use rokko
  use mod_frank
  implicit none

  integer::dim
  type(distributed_matrix)::mat,Z !defined in rokko
  type(grid)::g !defined in rokko
  type(solver)::solver_ !defined in rokko
  type(timer)::timer_
  
  real(8),allocatable::w(:),vec(:) !localized_vector
  character(len=100)::solver_name
  character(len=100)::tmp_str
  integer args_cnt, arg_len, status
  integer(c_int) :: n
  common /mydata/n

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
     solver_name = 'scalapack' !trim(tmp_str)
     call get_command_argument (2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  else
     write(*,*) "error: command line"
     write(*,*) "usage: solver_name matrix_size"
     stop
  endif

  call set_timer(timer_)
  call registrate_timer(timer_, 1, "diagonalize")

  write(*,*) "solver_name=", solver_name  
  write(*,*) "dim=",dim
!  dim = 10
  n = dim

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

  do count = 1, 1
    call generate_distributed_matrix_function(mat, c_funloc(func))

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    call diagonalize(solver_, mat, w, Z, timer_)
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  enddo
  write(*,*) "finised matrix generation frank"

  if (myrank_g .eq. 0) then
    write(*,*) "Computed Eigenvalues = "
    do i = 1, dim
      write(*,"(f30.20)") w(i)
    enddo
  endif

  if (myrank_g .eq. 0) then
    write(*,*) "num_procs = ", nprocs, nprocs_g
    write(*,*) "num_threads = ", omp_get_max_threads()
    write(*,*) "solver_name = ", trim(solver_name)
    write(*,*) "matrix = frank"
    write(*,*) "dim = ",dim
    write(*,*) "time = ", get_average_timer(timer_, 1)
  endif

  call del_distributed_matrix(mat)
  call del_distributed_matrix(Z)
  call del_timer(timer_)
  call del_solver(solver_)
  deallocate(w)

!!$  allocate(vec(dim))
!!$  do i=1, dim
!!$    call get_column_from_distributed_matrix(vec, Z, i)
!!$    write(*,"(a,i,a,100e25.15)") "Eigen vector of ",i, "=", vec
!!$  end do    
!!$  deallocate(vec)
  call MPI_finalize(ierr)
end program frank_matrix

