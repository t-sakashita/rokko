!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module frank_mod
  implicit none
  public frank_matrix_element
  integer, private :: dim

contains

  function frank_matrix_element(i, j)
    double precision :: frank_matrix_element
    integer, intent(in) :: i, j
    frank_matrix_element = dble(dim - max(i, j))
  end function frank_matrix_element

  subroutine frank_matrix_set_dimension(dim_in)
    integer, intent(in) :: dim_in
    dim = dim_in
  end subroutine frank_matrix_set_dimension

end module frank_mod

program frank_function
  use rokko
  use frank_mod
  use mpi
  implicit none
  integer :: dim
  type(rokko_parallel_dense_ev) :: solver
  type(rokko_grid) :: grid
  type(rokko_mapping_bc) :: map
  type(rokko_distributed_matrix) :: mat, Z
  type(rokko_eigen_vector) :: w
  character(len=:), allocatable :: library, routine
  character(len=:), allocatable :: library_routine, tmp_str

  integer :: provided,ierr, myrank, nprocs
  integer :: i

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() >= 1) then
     call get_command_argument_deferred(1, library_routine)
  else
     call rokko_parallel_dense_ev_default_solver(library_routine)
  endif
  call rokko_split_solver_name(library_routine, library, routine)

  if (command_argument_count() == 2) then
     call get_command_argument_deferred(2, tmp_str)
     read(tmp_str, *) dim
  else
     dim = 10
  endif

  if (myrank == 0) then
     print *,"library = ", library
     print *,"routine = ", routine
     print *,"dimension = ", dim
  endif

  call rokko_construct(solver, library)
  call rokko_construct(grid, MPI_COMM_WORLD, rokko_grid_row_major)
  call rokko_default_mapping(solver, dim, grid, map)
  call rokko_construct(mat, map)
  call rokko_construct(Z, map)
  call rokko_construct(w, dim)

  ! generate frank matrix from frank_matrix_element function
  call frank_matrix_set_dimension(dim)
  call rokko_generate0(mat, frank_matrix_element)
  call rokko_print(mat)

  call rokko_diagonalize(solver, mat, w, Z)

  if (myrank.eq.0) then
     write(*,*) "Computed Eigenvalues = "
     do i = 1, dim
        write(*,"(f30.20)") rokko_get_elem_f(w, i)
     enddo
  endif

  call rokko_destruct(mat)
  call rokko_destruct(Z)
  call rokko_destruct(w)
  call rokko_destruct(solver)
  call rokko_destruct(grid)

  call MPI_finalize(ierr)
end program frank_function

