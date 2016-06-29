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

module rokko_parallel_sparse_ev_mod
  use iso_c_binding
  use rokko_distributed_crs_matrix_mod, only : rokko_distributed_crs_matrix
  use rokko_distributed_mfree_mod, only : rokko_distributed_mfree
  use parameters
  use rokko_string
  implicit none
  
  type, bind(c) :: rokko_parallel_sparse_ev
     type(c_ptr) :: ptr
  end type rokko_parallel_sparse_ev

  interface
     subroutine rokko_parallel_sparse_ev_destruct(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
     end subroutine rokko_parallel_sparse_ev_destruct

     function rokko_parallel_sparse_ev_num_conv(solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_parallel_sparse_ev_num_conv
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
     end function rokko_parallel_sparse_ev_num_conv

     function rokko_parallel_sparse_ev_eigenvalue(solver, i) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       real(c_double) :: rokko_parallel_sparse_ev_eigenvalue
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: i
     end function rokko_parallel_sparse_ev_eigenvalue

     subroutine rokko_parallel_sparse_ev_eigenvector(solver, i, eig_vec) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
       integer(c_int), value, intent(in) :: i
       real(c_double), dimension(*), intent(inout) :: eig_vec
     end subroutine rokko_parallel_sparse_ev_eigenvector
  end interface

  interface
     type(c_ptr) function rokko_parallel_sparse_ev_default_solver_c() &
          bind(c,name='rokko_parallel_sparse_ev_default_solver')
       use iso_c_binding
       implicit none
     end function rokko_parallel_sparse_ev_default_solver_c

     integer(c_int) function rokko_parallel_sparse_ev_num_solvers_c() &
          bind(c,name='rokko_parallel_sparse_ev_num_solvers')
       use iso_c_binding
       implicit none
     end function rokko_parallel_sparse_ev_num_solvers_c

     type(c_ptr) function rokko_parallel_sparse_ev_solvers_c() &
          bind(c,name='rokko_parallel_sparse_ev_solvers')
       use iso_c_binding
       implicit none
     end function rokko_parallel_sparse_ev_solvers_c
  end interface
  
  interface rokko_parallel_sparse_ev_diagonalize
      
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(solver, mat, params, params_out) &
          bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix
     
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_noreturn(solver, &
          mat, params) bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_noreturn_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_crs_matrix), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_crs_noreturn
     
     
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree(solver, mat, params, params_out) &
          bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_mfree_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_mfree
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_mfree), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
       type(rokko_parameters), intent(out) :: params_out
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree
     
     subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn(solver, mat, params) &
          bind(c,name="rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn_f")
       use iso_c_binding
       use parameters
       import rokko_parallel_sparse_ev, rokko_distributed_mfree
       implicit none
       type(rokko_parallel_sparse_ev), intent(inout) :: solver
       type(rokko_distributed_mfree), intent(inout) :: mat
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn

  end interface rokko_parallel_sparse_ev_diagonalize

  ! constructor for CRS
  interface
     subroutine rokko_distributed_crs_matrix_construct(matrix, dim1, dim2, solver) bind(c)
       use iso_c_binding
       import rokko_parallel_sparse_ev, rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       type(rokko_parallel_sparse_ev), value, intent(in) :: solver
     end subroutine rokko_distributed_crs_matrix_construct
  end interface

  ! generic names
  interface rokko_construct
     procedure rokko_parallel_sparse_ev_construct
     procedure rokko_distributed_crs_matrix_construct
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_parallel_sparse_ev_destruct
  end interface rokko_destruct

  interface rokko_diagonalize
     procedure rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix
     procedure rokko_parallel_sparse_ev_diagonalize_distributed_crs_noreturn
     procedure rokko_parallel_sparse_ev_diagonalize_distributed_mfree
     procedure rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn
  end interface rokko_diagonalize

  interface rokko_num_conv
     procedure rokko_parallel_sparse_ev_num_conv
  end interface rokko_num_conv

  interface rokko_eigenvalue
     procedure rokko_parallel_sparse_ev_eigenvalue
  end interface rokko_eigenvalue
  
  interface rokko_eigenvector
     procedure rokko_parallel_sparse_ev_eigenvector
  end interface rokko_eigenvector
  
contains

  subroutine rokko_parallel_sparse_ev_default_solver(name)
    character(len=*), intent(out) :: name
    type(c_ptr) :: name_ptr
    name_ptr = rokko_parallel_sparse_ev_default_solver_c ()
    call rokko_get_string_fixedsize (name_ptr, name)
  end subroutine rokko_parallel_sparse_ev_default_solver

  subroutine rokko_parallel_sparse_ev_num_solvers(num)
    integer, intent(out) :: num
    num = rokko_parallel_sparse_ev_num_solvers_c()
  end subroutine rokko_parallel_sparse_ev_num_solvers

  subroutine rokko_parallel_sparse_ev_solvers(names)
    implicit none
    type(string), allocatable, intent(out) :: names(:)
    type(c_ptr) :: ptr, ptr_i
    integer :: i, size
    character(len=:), allocatable :: str
    ptr = rokko_parallel_sparse_ev_solvers_c ()
    size = rokko_parallel_sparse_ev_num_solvers_c ()
    allocate(names(size))
    do i = 1, size
       ptr_i = rokko_string_i_c (ptr, i-1)
       call rokko_get_string(ptr_i, str)
       names(i)%str = str
    enddo
  end subroutine rokko_parallel_sparse_ev_solvers

  subroutine rokko_parallel_sparse_ev_construct(solver, solver_name)
    implicit none
    interface
       subroutine rokko_parallel_sparse_ev_construct_f(solver, solver_name) bind(c)
         use iso_c_binding
         import rokko_parallel_sparse_ev
         implicit none
         type(rokko_parallel_sparse_ev), intent(out) :: solver
         character(kind=c_char), intent(in) :: solver_name(*)
       end subroutine rokko_parallel_sparse_ev_construct_f
    end interface
    type(rokko_parallel_sparse_ev), intent(inout) :: solver
    character(*), intent(in) :: solver_name
    call rokko_parallel_sparse_ev_construct_f(solver, trim(solver_name)//C_NULL_CHAR)
  end subroutine rokko_parallel_sparse_ev_construct

end module rokko_parallel_sparse_ev_mod

