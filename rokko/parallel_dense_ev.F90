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

#include <rokko/config.h>

module rokko_parallel_dense_ev_mod
  use iso_c_binding
  use rokko_grid_mod
  use rokko_mapping_bc_mod, only : rokko_mapping_bc
  use rokko_distributed_matrix_mod, only : rokko_distributed_matrix
  use rokko_serial_dense
  use rokko_string
  implicit none

  type, bind(c) :: rokko_parallel_dense_ev
     type(c_ptr) ptr
  end type rokko_parallel_dense_ev
  
  interface
     subroutine rokko_parallel_dense_ev_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_dense_ev
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
     end subroutine rokko_parallel_dense_ev_destruct

     type(c_ptr) function rokko_parallel_dense_ev_default_solver_c() &
          bind(c,name='rokko_parallel_dense_ev_default_solver')
       use iso_c_binding
       implicit none
     end function rokko_parallel_dense_ev_default_solver_c

     integer(c_int) function rokko_parallel_dense_ev_num_solvers_c() &
          bind(c,name='rokko_parallel_dense_ev_num_solvers')
       use iso_c_binding
       implicit none
     end function rokko_parallel_dense_ev_num_solvers_c

     type(c_ptr) function rokko_parallel_dense_ev_solvers_c() &
          bind(c,name='rokko_parallel_dense_ev_solvers')
       use iso_c_binding
       implicit none
     end function rokko_parallel_dense_ev_solvers_c
  end interface

  interface
     subroutine rokko_parallel_dense_ev_default_mapping(solver, global_dim, grid, map) &
          bind(c,name="rokko_parallel_dense_ev_default_mapping_f")
       use iso_c_binding
       import rokko_parallel_dense_ev, rokko_grid, rokko_mapping_bc
       implicit none
       type(rokko_parallel_dense_ev), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: global_dim
       type(rokko_grid), value, intent(in) :: grid
       type(rokko_mapping_bc), intent(out) :: map
     end subroutine rokko_parallel_dense_ev_default_mapping
  end interface
  
  interface rokko_parallel_dense_ev_diagonalize

     subroutine rokko_parallel_dense_ev_diagonalize(solver, mat, &
          eigvals, eigvecs, params, params_out) bind(c,name="rokko_parallel_dense_ev_diagonalize_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_dense_ev, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_matrix), intent(inout) :: eigvecs
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_parallel_dense_ev_diagonalize
     
     subroutine rokko_parallel_dense_ev_diagonalize_no_params_out(solver, mat, &
          eigvals, eigvecs, params) bind(c,name="rokko_parallel_dense_ev_diagonalize_no_params_out_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_dense_ev, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_matrix), intent(inout) :: eigvecs
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parallel_dense_ev_diagonalize_no_params_out
          
     subroutine rokko_parallel_dense_ev_diagonalize_no_params_inout(solver, mat, &
          eigvals, eigvecs) bind(c,name='rokko_parallel_dense_ev_diagonalize_no_params_inout_f')
       use iso_c_binding
       use parameters
       import rokko_parallel_dense_ev, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_matrix), intent(inout) :: eigvecs
     end subroutine rokko_parallel_dense_ev_diagonalize_no_params_inout

     ! Only eigenvalues
     subroutine rokko_parallel_dense_ev_diagonalize_eigvals(solver, mat, eigvals, params, params_out) &
          bind(c,name="rokko_parallel_dense_ev_diagonalize_eigvals_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_dense_ev, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_parallel_dense_ev_diagonalize_eigvals
     
     subroutine rokko_parallel_dense_ev_diagonalize_eigvals_no_params_out(solver, mat, eigvals, params) &
          bind(c,name="rokko_parallel_dense_ev_diagonalize_eigvals_no_params_out_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_dense_ev, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parallel_dense_ev_diagonalize_eigvals_no_params_out
     
     subroutine rokko_parallel_dense_ev_diagonalize_eigvals_no_params_inout(solver, mat, eigvals) &
          bind(c,name="rokko_parallel_dense_ev_diagonalize_eigvals_no_params_inout_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_dense_ev, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_parallel_dense_ev), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
     end subroutine rokko_parallel_dense_ev_diagonalize_eigvals_no_params_inout     
     
  end interface rokko_parallel_dense_ev_diagonalize

  ! generic names
  !interface construct
  !   procedure rokko_parallel_dense_ev_construct
  !end interface contruct

  interface destruct
     procedure rokko_parallel_dense_ev_destruct
  end interface destruct

  interface diagonalize
     procedure rokko_parallel_dense_ev_diagonalize
  end interface diagonalize

  !
  ! rokko_frank_matrix for parallel dense solvers
  !

  ! generic name
  interface frank_matrix_generate
     procedure rokko_frank_matrix_generate_distributed_matrix
  end interface frank_matrix_generate
  
  interface
     subroutine rokko_frank_matrix_generate_distributed_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end subroutine rokko_frank_matrix_generate_distributed_matrix
  end interface

contains

  subroutine rokko_parallel_dense_ev_default_solver(name)
    use rokko_string
    character(len=*), intent(out) :: name
    type(c_ptr) :: name_ptr
    name_ptr = rokko_parallel_dense_ev_default_solver_c ()
    call rokko_get_string_fixedsize (name_ptr, name)
  end subroutine rokko_parallel_dense_ev_default_solver

  subroutine rokko_parallel_dense_ev_num_solvers(num)
    integer, intent(out) :: num
    num = rokko_parallel_dense_ev_num_solvers_c()
  end subroutine rokko_parallel_dense_ev_num_solvers

  subroutine rokko_parallel_dense_ev_solvers(names)
    use iso_c_binding
    implicit none
    type(string), allocatable, intent(out) :: names(:)
    type(c_ptr) :: ptr, ptr_i
    integer :: i, size
    character(len=:), allocatable :: str
    ptr = rokko_parallel_dense_ev_solvers_c ()
    size = rokko_parallel_dense_ev_num_solvers_c ()
    allocate(names(size))
    do i = 1, size
       ptr_i = rokko_string_i_c (ptr, i-1)
       call rokko_get_string(ptr_i, str)
       names(i)%str = str
    enddo
  end subroutine rokko_parallel_dense_ev_solvers
  
  subroutine rokko_parallel_dense_ev_construct(solver, solver_name)
    use iso_c_binding
    use rokko_parallel_dense_ev_type
    implicit none
    interface
       subroutine rokko_parallel_dense_ev_construct_f(solver, solver_name) bind(c)
         use rokko_parallel_dense_ev_type
         implicit none
         type(rokko_parallel_dense_ev), intent(out) :: solver
         character(kind=c_char), intent(in) :: solver_name(*)
       end subroutine rokko_parallel_dense_ev_construct_f
    end interface
    type(rokko_parallel_dense_ev), intent(inout) :: solver
    character(*), intent(in) :: solver_name
    call rokko_parallel_dense_ev_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
  end subroutine rokko_parallel_dense_ev_construct
  
end module rokko_parallel_dense_ev_mod

