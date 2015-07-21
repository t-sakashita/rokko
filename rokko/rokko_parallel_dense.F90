!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

#include <rokko/config.h>

module rokko_parallel_dense
  use iso_c_binding
  use rokko_parallel_dense_classes
  use rokko_distributed_matrix_mod
  implicit none
  
  !
  ! rokko_grid
  !
  
  interface
     subroutine rokko_grid_destruct(grid) bind(c)
       import rokko_grid
       implicit none
       type(rokko_grid), intent(inout) :: grid
     end subroutine rokko_grid_destruct

     function rokko_grid_get_myrank(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       integer(c_int) :: rokko_grid_get_myrank
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_myrank

     function rokko_grid_get_nprocs(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       integer(c_int) :: rokko_grid_get_nprocs
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_nprocs
  end interface
  
  interface rokko_grid_construct
     subroutine rokko_grid_construct_f(grid, comm, grid_major) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       type(rokko_grid), intent(out) :: grid
       integer(c_int), value, intent(in) :: comm
       integer(c_int), value, intent(in) :: grid_major
     end subroutine rokko_grid_construct_f
  end interface rokko_grid_construct

  !
  ! rokko_parallel_dense_solver
  !
     
  interface
     subroutine rokko_parallel_dense_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_dense_solver
       implicit none
       type(rokko_parallel_dense_solver), intent(inout) :: solver
     end subroutine rokko_parallel_dense_solver_destruct
     
     subroutine rokko_parallel_dense_solver_diagonalize_distributed_matrix(solver, mat, eigvals, &
          eigvecs) bind(c)
       use iso_c_binding
       use rokko_serial_dense, only : rokko_localized_vector
       import rokko_parallel_dense_solver, rokko_distributed_matrix
       implicit none
       type(rokko_parallel_dense_solver), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_matrix), intent(inout) :: eigvecs
     end subroutine rokko_parallel_dense_solver_diagonalize_distributed_matrix
  end interface
  
  !
  ! collective operations
  !

  interface
     function rokko_gather(matrix, array, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_gather
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value ::root
     end function rokko_gather

     subroutine rokko_all_gather(matrix, array)
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       double precision, intent(in), target :: array(:,:)
     end subroutine rokko_all_gather

     function rokko_scatter(array, matrix, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_scatter
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value :: root
     end function rokko_scatter
  end interface

  !
  ! rokko_frank_matrix for parallel dense solvers
  !

  interface
     subroutine rokko_frank_matrix_generate_distributed_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_frank_matrix_generate_distributed_matrix
  end interface
  
end module rokko_parallel_dense
