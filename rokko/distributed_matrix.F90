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

module rokko_distributed_matrix_mod
  use iso_c_binding
  use rokko_grid_mod, only : rokko_grid
  use rokko_mapping_bc_mod, only : rokko_mapping_bc
  implicit none

  type, bind(c) :: rokko_distributed_matrix
     type(c_ptr) :: ptr
     integer(c_int) :: major
  end type rokko_distributed_matrix

  ! generic names
  interface rokko_construct
     procedure rokko_distributed_matrix_construct
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_distributed_matrix_destruct
  end interface rokko_destruct

  interface rokko_get_global
     procedure rokko_distributed_matrix_get_global
  end interface rokko_get_global

  interface rokko_get_local
     procedure rokko_distributed_matrix_get_local
  end interface rokko_get_local
  
  interface rokko_set_global
     procedure rokko_distributed_matrix_set_global
  end interface rokko_set_global

  interface rokko_set_local
     procedure rokko_distributed_matrix_set_local
  end interface rokko_set_local

  interface rokko_get_m_global
     procedure rokko_distributed_matrix_get_m_global
  end interface rokko_get_m_global

  interface rokko_get_n_global
     procedure rokko_distributed_matrix_get_n_global
  end interface rokko_get_n_global
   
  interface rokko_get_m_local
     procedure rokko_distributed_matrix_get_m_local
  end interface rokko_get_m_local

  interface rokko_get_n_local
     procedure rokko_distributed_matrix_get_n_local
  end interface rokko_get_n_local

  interface rokko_translate_l2g_row
     procedure  rokko_distributed_matrix_translate_l2g_row
  end interface rokko_translate_l2g_row

  interface rokko_translate_l2g_col
     procedure  rokko_distributed_matrix_translate_l2g_col
  end interface rokko_translate_l2g_col

  interface rokko_translate_g2l_row
     procedure  rokko_distributed_matrix_translate_g2l_row
  end interface rokko_translate_g2l_row

  interface rokko_translate_g2l_col
     procedure  rokko_distributed_matrix_translate_g2l_col
  end interface rokko_translate_g2l_col

  interface rokko_generate
     module procedure rokko_distributed_matrix_generate_function
     module procedure rokko_distributed_matrix_generate_from_array
  end interface rokko_generate

  interface rokko_get_array_pointer
     module procedure rokko_distributed_matrix_get_array_pointer
  end interface rokko_get_array_pointer

  interface rokko_print
     procedure rokko_distributed_matrix_print
  end interface rokko_print

  interface
     subroutine rokko_distributed_matrix_construct(matrix, map) &
          bind(c)
       use iso_c_binding
       import rokko_mapping_bc, rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(rokko_mapping_bc), value, intent(in) :: map
     end subroutine rokko_distributed_matrix_construct
     
     subroutine rokko_distributed_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_matrix_destruct
     
     subroutine rokko_distributed_matrix_generate_function_c(matrix, cproc) &
          bind(c,name='rokko_distributed_matrix_generate_function')
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_distributed_matrix_generate_function_c
     
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
       type(rokko_distributed_matrix), value, intent(in) :: matrix       
       integer(c_int), value, intent(in) :: local_i, local_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_local
     
     function rokko_distributed_matrix_get_local(matrix, local_i,local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i,local_j
     end function rokko_distributed_matrix_get_local
     
     subroutine rokko_distributed_matrix_set_global(matrix, global_i, global_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_i, global_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_global
     
     function rokko_distributed_matrix_get_global(matrix, global_i, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in):: global_i, global_j
     end function rokko_distributed_matrix_get_global
     
     function rokko_distributed_matrix_get_m_local(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_m_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_local

     function rokko_distributed_matrix_get_n_local(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_n_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_local

     function rokko_distributed_matrix_get_m_global(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_m_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_global

     function rokko_distributed_matrix_get_n_global(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_n_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_global

     function rokko_distributed_matrix_get_nprocs(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_nprocs
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_nprocs

     function rokko_distributed_matrix_get_myrank(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_myrank
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_myrank

     function rokko_distributed_matrix_get_nprow(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_nprow
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_nprow

     function rokko_distributed_matrix_get_npcol(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_npcol
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_npcol

     function rokko_distributed_matrix_get_myrow(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_myrow
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_myrow

     function rokko_distributed_matrix_get_mycol(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_mycol
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_mycol

     function rokko_distributed_matrix_translate_l2g_row(matrix, local_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_row
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value,intent(in) :: local_i
     end function rokko_distributed_matrix_translate_l2g_row

     function rokko_distributed_matrix_translate_l2g_col(matrix, local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_col
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::local_j
     end function rokko_distributed_matrix_translate_l2g_col

     function rokko_distributed_matrix_translate_g2l_row(matrix, global_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_row
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::global_i
     end function rokko_distributed_matrix_translate_g2l_row

     function rokko_distributed_matrix_translate_g2l_col(matrix, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_col
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_j
     end function rokko_distributed_matrix_translate_g2l_col

     function rokko_distributed_matrix_is_row_major(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       logical(c_bool) :: rokko_distributed_matrix_is_row_major
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_is_row_major

     function rokko_distributed_matrix_is_col_major(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       logical(c_bool) :: rokko_distributed_matrix_is_col_major
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_is_col_major
     
     type(c_ptr) function rokko_distributed_matrix_get_array_pointer_c(matrix) &
          bind(c,name='rokko_distributed_matrix_get_array_pointer')
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_ptr) :: c_array_ptr
     end function rokko_distributed_matrix_get_array_pointer_c
  end interface

contains
  
  subroutine rokko_distributed_matrix_generate_from_array(matrix, array)
    type(rokko_distributed_matrix), intent(out) :: matrix
    double precision, intent(in) :: array(:,:)
    integer :: m_local, n_local, local_i, local_j, global_i, global_j
    m_local = rokko_distributed_matrix_get_m_local(matrix)
    n_local = rokko_distributed_matrix_get_n_local(matrix)
    do local_i = 0, m_local-1
       do local_j = 0, n_local-1
          global_i = rokko_distributed_matrix_translate_l2g_row(matrix, local_i)
          global_j = rokko_distributed_matrix_translate_l2g_col(matrix, local_j)
          call rokko_distributed_matrix_set_local(matrix, local_i, local_j, &
               array(global_i+1, global_j+1))
       enddo
    enddo
  end subroutine rokko_distributed_matrix_generate_from_array

  subroutine rokko_distributed_matrix_generate_function(matrix, func_in)
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    type(c_funptr) :: cproc
    interface
       double precision function func_in (i, j) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int), value, intent(in) :: i, j
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_distributed_matrix_generate_function_c(matrix, cproc)
  end subroutine rokko_distributed_matrix_generate_function

  subroutine rokko_distributed_matrix_get_array_pointer(matrix, f_array_ptr)
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    double precision, pointer, dimension(:,:), intent(out) :: f_array_ptr
    type(c_ptr) :: c_array_ptr
    integer(c_int) :: m_local, n_local
    c_array_ptr = rokko_distributed_matrix_get_array_pointer_c(matrix)
    m_local = rokko_distributed_matrix_get_m_local(matrix)
    n_local = rokko_distributed_matrix_get_n_local(matrix)
    call c_f_pointer(c_array_ptr, f_array_ptr, (/m_local,n_local/) )
  end subroutine rokko_distributed_matrix_get_array_pointer
  
end module rokko_distributed_matrix_mod

