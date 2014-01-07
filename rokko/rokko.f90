!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>,
!                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module rokko
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_grid_col_major = 1, rokko_grid_row_major = 2
     enumerator :: rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
  end enum

  type, bind(c) :: rokko_grid
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_grid

  type, bind(c) :: rokko_solver
     type(c_ptr) ptr
  end type rokko_solver
  
  type, bind(c) :: rokko_localized_vector
     type(c_ptr) ptr
  end type rokko_localized_vector

  type, bind(c) :: rokko_localized_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_localized_matrix

  type, bind(c) :: rokko_distributed_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_distributed_matrix

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
  
  interface

     ! rokko_grid
     
     subroutine rokko_grid_destruct(grid) bind(c)
       import rokko_grid
       implicit none
       type(rokko_grid), intent(inout) :: grid
     end subroutine rokko_grid_destruct

     integer(c_int) function rokko_grid_get_myrank(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_myrank

     integer(c_int) function rokko_grid_get_nprocs(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_nprocs

     ! rokko_solver

     subroutine rokko_solver_construct(solver, solver_name)
       import rokko_solver
       implicit none
       type(rokko_solver), intent(inout) :: solver
       character(*), intent(in) :: solver_name
     end subroutine rokko_solver_construct
     
     subroutine rokko_solver_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_solver
       implicit none
       type(rokko_solver), intent(inout) :: solver
     end subroutine rokko_solver_destruct
     
     subroutine rokko_solver_diagonalize_distributed_matrix(solver, mat, eigvals, eigvecs) bind(c)
       use iso_c_binding
       import rokko_solver, rokko_distributed_matrix, rokko_localized_vector
       implicit none
       type(rokko_solver), intent(inout) :: solver
       type(rokko_distributed_matrix), intent(inout) :: mat
       type(rokko_localized_vector), intent(inout) :: eigvals
       type(rokko_distributed_matrix), intent(inout) :: eigvecs
     end subroutine rokko_solver_diagonalize_distributed_matrix

     ! rokko_localized_vector

     subroutine rokko_localized_vector_construct(vec, dim1) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), intent(out) :: vec
       integer(c_int), value, intent(in) :: dim1
     end subroutine rokko_localized_vector_construct

     subroutine rokko_localized_vector_destruct(vec) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), intent(inout) :: vec
     end subroutine rokko_localized_vector_destruct
  end interface

  interface rokko_localized_vector_get
     real(c_double) function rokko_localized_vector_get_f(vec, i) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
     end function rokko_localized_vector_get_f
  end interface rokko_localized_vector_get

  interface
     
     ! rokko_localized_matrix

     subroutine rokko_localized_matrix_construct(matrix, dim1, dim2, matrix_major) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       integer(c_int), value, intent(in) :: matrix_major
     end subroutine rokko_localized_matrix_construct
     
     subroutine rokko_localized_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), intent(inout) :: matrix
     end subroutine rokko_localized_matrix_destruct
     
     ! rokko_distributed_matrix

     subroutine rokko_distributed_matrix_construct(matrix, dim1, dim2, grid, solver, matrix_major) &
          bind(c)
       use iso_c_binding
       import rokko_grid, rokko_solver, rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       type(rokko_grid), value, intent(in) :: grid
       type(rokko_solver), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: matrix_major
     end subroutine rokko_distributed_matrix_construct
     
     subroutine rokko_distributed_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_matrix_destruct

     subroutine rokko_distributed_matrix_generate_function(matrix, func) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(c_funptr), intent(in), value :: func
     end subroutine rokko_distributed_matrix_generate_function
     
     subroutine rokko_distributed_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_matrix_print
     
     subroutine rokko_distributed_matrix_set_local(matrix, local_i, local_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int),value, intent(in):: local_i,local_j
       real(c_double),value, intent(in):: value
     end subroutine rokko_distributed_matrix_set_local

     real(c_double) function rokko_distributed_matrix_get_local(matrix, local_i,local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value, intent(in):: local_i,local_j
     end function rokko_distributed_matrix_get_local

     subroutine rokko_distributed_matrix_set_global(matrix, global_i,&
        & global_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int),value, intent(in):: global_i,global_j
       real(c_double),value, intent(in):: value
     end subroutine rokko_distributed_matrix_set_global
     
     real(c_double) function rokko_distributed_matrix_get_global(matrix, global_i,global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value, intent(in):: global_i,global_j
     end function rokko_distributed_matrix_get_global

     integer(c_int) function rokko_distributed_matrix_get_m_local(matrix)&
        & bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_local

     integer(c_int) function rokko_distributed_matrix_get_n_local(matrix)&
        & bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_local

     integer(c_int) function rokko_distributed_matrix_get_m_global(matrix)&
        & bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_global

     integer(c_int) function rokko_distributed_matrix_get_n_global(matrix)&
        & bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_global

     integer(c_int) function rokko_distributed_matrix_get_nprocs(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_nprocs

     integer(c_int) function rokko_distributed_matrix_get_myrank(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_myrank
       

     integer(c_int) function rokko_distributed_matrix_translate_l2g_row(matrix,local_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
       integer(c_int),value,intent(in)::local_i
     end function rokko_distributed_matrix_translate_l2g_row

     integer(c_int) function rokko_distributed_matrix_translate_l2g_col(matrix, local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
       integer(c_int),value,intent(in)::local_j
     end function rokko_distributed_matrix_translate_l2g_col


     integer(c_int) function rokko_distributed_matrix_translate_g2l_row(matrix, global_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
       integer(c_int),value,intent(in)::global_i
     end function rokko_distributed_matrix_translate_g2l_row

     integer(c_int) function rokko_distributed_matrix_translate_g2l_col(matrix, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),value, intent(in) :: matrix
       integer(c_int),value,intent(in)::global_j
     end function rokko_distributed_matrix_translate_g2l_col


!!$     subroutine rokko_distributed_matrix_generate_array(matrix,&
!!$        & array) 
!!$       import rokko_distributed_matrix
!!$       implicit none
!!$       type(rokko_distributed_matrix),intent(out) :: matrix
!!$       real*8,intent(in):: array(:,:)
!!$     end subroutine rokko_distributed_matrix_generate_array


     ! rokko_cellective
     integer(c_int) function rokko_gather(matrix, array, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),intent(out) ::matrix
       type(c_ptr),value,intent(in)::array
       integer(c_int),value::root
     end function rokko_gather

     integer(c_int) function rokko_scatter(array, matrix, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix),intent(out) ::matrix
       type(c_ptr),value,intent(in)::array
       integer(c_int),value::root
     end function rokko_scatter
       
  end interface
contains
  subroutine rokko_distributed_matrix_generate_array(matrix, array)
    implicit none
    type(rokko_distributed_matrix),intent(out) :: matrix
    real*8,intent(in):: array(:,:)
    integer:: m_local,n_local,local_i,local_j,global_i,global_j
    m_local = rokko_distributed_matrix_get_m_local(matrix)
    n_local = rokko_distributed_matrix_get_n_local(matrix)
    do local_i = 0, m_local - 1 
      do local_j = 0, n_local -1
        global_i = rokko_distributed_matrix_translate_l2g_row(matrix, local_i)
        global_j = rokko_distributed_matrix_translate_l2g_col(matrix, local_j)

        call rokko_distributed_matrix_set_local(matrix, local_i,&
           & local_j, array(global_i + 1, global_j + 1))
      enddo
    enddo
  end subroutine rokko_distributed_matrix_generate_array
  subroutine rokko_all_gather(matrix, array)
    implicit none
    type(rokko_distributed_matrix),intent(out) ::matrix
    real*8,intent(in),target:: array(:,:)
    integer(c_int)::root, nprocs,ierr
    nprocs = rokko_distributed_matrix_get_nprocs(matrix)

    do root=0, nprocs - 1
      ierr = rokko_gather(matrix, c_loc(array(1,1)), root)
    end do
  end subroutine rokko_all_gather

end module rokko
