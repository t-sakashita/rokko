!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
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
     module procedure rokko_distributed_matrix_construct_array1
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_distributed_matrix_destruct
  end interface rokko_destruct

  interface rokko_get_global0
     procedure rokko_distributed_matrix_get_global0
  end interface rokko_get_global0

  interface rokko_get_local0
     procedure rokko_distributed_matrix_get_local0
  end interface rokko_get_local0

  interface rokko_set_global0
     procedure rokko_distributed_matrix_set_global0
  end interface rokko_set_global0

  interface rokko_set_local0
     procedure rokko_distributed_matrix_set_local0
  end interface rokko_set_local0

  ! offset by one
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

  interface rokko_get_mb
     procedure rokko_distributed_matrix_get_mb
  end interface rokko_get_mb

  interface rokko_get_nb
     procedure rokko_distributed_matrix_get_nb
  end interface rokko_get_nb

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

  interface rokko_get_lld
     procedure rokko_distributed_matrix_get_lld
  end interface rokko_get_lld

  interface rokko_get_m_size
     procedure rokko_distributed_matrix_get_m_size
  end interface rokko_get_m_size

  interface rokko_get_n_size
     procedure rokko_distributed_matrix_get_n_size
  end interface rokko_get_n_size

  interface rokko_translate_l2g_row0
     procedure rokko_distributed_matrix_translate_l2g_row0
  end interface rokko_translate_l2g_row0

  interface rokko_translate_l2g_col0
     procedure rokko_distributed_matrix_translate_l2g_col0
  end interface rokko_translate_l2g_col0

  interface rokko_translate_g2l_row0
     procedure rokko_distributed_matrix_translate_g2l_row0
  end interface rokko_translate_g2l_row0

  interface rokko_translate_g2l_col0
     procedure rokko_distributed_matrix_translate_g2l_col0
  end interface rokko_translate_g2l_col0

  ! offset by one
  interface rokko_translate_l2g_row
     procedure rokko_distributed_matrix_translate_l2g_row
  end interface rokko_translate_l2g_row

  interface rokko_translate_l2g_col
     procedure rokko_distributed_matrix_translate_l2g_col
  end interface rokko_translate_l2g_col

  interface rokko_translate_g2l_row
     procedure rokko_distributed_matrix_translate_g2l_row
  end interface rokko_translate_g2l_row

  interface rokko_translate_g2l_col
     procedure rokko_distributed_matrix_translate_g2l_col
  end interface rokko_translate_g2l_col

  interface rokko_is_row_major
     procedure rokko_distributed_matrix_is_row_major
  end interface rokko_is_row_major

  interface rokko_is_col_major
     procedure rokko_distributed_matrix_is_col_major
  end interface rokko_is_col_major

  interface rokko_generate0
     module procedure rokko_distributed_matrix_generate_function0
  end interface rokko_generate0

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

     subroutine rokko_distributed_matrix_construct_array_sizes(matrix, map, dim1, dim2, array) &
          bind(c)
       use iso_c_binding
       import rokko_mapping_bc, rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in) :: dim1, dim2
       double precision, intent(in) :: array(dim1, dim2)
     end subroutine rokko_distributed_matrix_construct_array_sizes

     subroutine rokko_distributed_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_matrix_destruct
     
     subroutine rokko_distributed_matrix_generate_function0_p(matrix, cproc) &
          & bind(c,name="rokko_distributed_matrix_generate_function_p")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_distributed_matrix_generate_function0_p

     subroutine rokko_distributed_matrix_generate_function1_p(matrix, cproc) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_distributed_matrix_generate_function1_p

     subroutine rokko_distributed_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_matrix_print

     subroutine rokko_distributed_matrix_set_local0(matrix, local_i, local_j, value) &
          & bind(c,name="rokko_distributed_matrix_set_local")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i, local_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_local0

     function rokko_distributed_matrix_get_local0(matrix, local_i,local_j) &
          & bind(c,name="rokko_distributed_matrix_get_local")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_local0
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i,local_j
     end function rokko_distributed_matrix_get_local0

     subroutine rokko_distributed_matrix_set_global0(matrix, global_i, global_j, value) &
          & bind(c,name="rokko_distributed_matrix_set_global")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_i, global_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_global0
     
     function rokko_distributed_matrix_get_global0(matrix, global_i, global_j) &
          & bind(c,name="rokko_distributed_matrix_get_global")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_global0
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in):: global_i, global_j
     end function rokko_distributed_matrix_get_global0

     ! offset by one
     subroutine rokko_distributed_matrix_set_local(matrix, local_i, local_j, value) &
          & bind(c,name="rokko_distributed_matrix_set_local1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i, local_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_local

     function rokko_distributed_matrix_get_local(matrix, local_i,local_j) &
          & bind(c,name="rokko_distributed_matrix_get_local1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i,local_j
     end function rokko_distributed_matrix_get_local

     subroutine rokko_distributed_matrix_set_global(matrix, global_i, global_j, value) &
          & bind(c,name="rokko_distributed_matrix_set_global1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_i, global_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_global

     function rokko_distributed_matrix_get_global(matrix, global_i, global_j) &
          & bind(c,name="rokko_distributed_matrix_get_global1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in):: global_i, global_j
     end function rokko_distributed_matrix_get_global

     function rokko_distributed_matrix_get_mb(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_mb
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_mb

     function rokko_distributed_matrix_get_nb(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_nb
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_nb

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

     function rokko_distributed_matrix_get_lld(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_lld
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_lld

     function rokko_distributed_matrix_get_m_size(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_m_size
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_size

     function rokko_distributed_matrix_get_n_size(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_n_size
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_size

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

     function rokko_distributed_matrix_translate_l2g_row0(matrix, local_i) &
       & bind(c,name="rokko_distributed_matrix_translate_l2g_row")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_row0
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value,intent(in) :: local_i
     end function rokko_distributed_matrix_translate_l2g_row0

     function rokko_distributed_matrix_translate_l2g_col0(matrix, local_j) &
          & bind(c,name="rokko_distributed_matrix_translate_l2g_col")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_col0
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::local_j
     end function rokko_distributed_matrix_translate_l2g_col0

     function rokko_distributed_matrix_translate_g2l_row0(matrix, global_i) &
          & bind(c,name="rokko_distributed_matrix_translate_g2l_row")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_row0
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::global_i
     end function rokko_distributed_matrix_translate_g2l_row0

     function rokko_distributed_matrix_translate_g2l_col0(matrix, global_j) &
          & bind(c,name="rokko_distributed_matrix_translate_g2l_col")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_col0
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_j
     end function rokko_distributed_matrix_translate_g2l_col0

     ! offset by one
     function rokko_distributed_matrix_translate_l2g_row(matrix, local_i) &
          & bind(c,name="rokko_distributed_matrix_translate_l2g_row1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_row
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value,intent(in) :: local_i
     end function rokko_distributed_matrix_translate_l2g_row

     function rokko_distributed_matrix_translate_l2g_col(matrix, local_j) &
          & bind(c,name="rokko_distributed_matrix_translate_l2g_col1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_col
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::local_j
     end function rokko_distributed_matrix_translate_l2g_col

     function rokko_distributed_matrix_translate_g2l_row(matrix, global_i) &
          & bind(c,name="rokko_distributed_matrix_translate_g2l_row1")
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_row
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::global_i
     end function rokko_distributed_matrix_translate_g2l_row

     function rokko_distributed_matrix_translate_g2l_col(matrix, global_j) &
          & bind(c,name="rokko_distributed_matrix_translate_g2l_col1")
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
          & bind(c,name='rokko_distributed_matrix_get_array_pointer')
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_ptr) :: c_array_ptr
     end function rokko_distributed_matrix_get_array_pointer_c
  end interface

contains

  subroutine rokko_distributed_matrix_construct_array1(matrix, map, array) bind(c)
    use iso_c_binding
    implicit none
    type(rokko_distributed_matrix), intent(out) :: matrix
    type(rokko_mapping_bc), value, intent(in) :: map
    double precision, intent(in) :: array(:,:)
    integer :: sizes(2)

    sizes = shape(array)
    call rokko_distributed_matrix_construct_array_sizes(matrix, map, sizes(1), sizes(2), array)
  end subroutine rokko_distributed_matrix_construct_array1

  subroutine rokko_distributed_matrix_generate_from_array(matrix, array)
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    double precision, intent(in) :: array(:,:)
    integer :: m_local, n_local, local_i, local_j, global_i, global_j
    m_local = rokko_distributed_matrix_get_m_local(matrix)
    n_local = rokko_distributed_matrix_get_n_local(matrix)
    do local_i = 0, m_local-1
       global_i = rokko_distributed_matrix_translate_l2g_row0(matrix, local_i)
       do local_j = 0, n_local-1
          global_j = rokko_distributed_matrix_translate_l2g_col0(matrix, local_j)
          call rokko_distributed_matrix_set_local0(matrix, local_i, local_j, &
               & array(global_i+1, global_j+1))
       enddo
    enddo
  end subroutine rokko_distributed_matrix_generate_from_array

  subroutine rokko_distributed_matrix_generate_function0(matrix, func_in)
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    type(c_funptr) :: cproc
    interface
       function func_in (i, j)
         use, intrinsic :: iso_c_binding
         double precision :: func_in
         integer, intent(in) :: i, j
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_distributed_matrix_generate_function0_p(matrix, cproc)
  end subroutine rokko_distributed_matrix_generate_function0

  subroutine rokko_distributed_matrix_generate_function(matrix, func_in)
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    type(c_funptr) :: cproc
    interface
       function func_in (i, j)
         implicit none
         double precision :: func_in
         integer, intent(in) :: i, j
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_distributed_matrix_generate_function1_p(matrix, cproc)
  end subroutine rokko_distributed_matrix_generate_function

  subroutine rokko_distributed_matrix_get_array_pointer(matrix, f_array_ptr)
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    double precision, pointer, dimension(:,:), intent(out) :: f_array_ptr
    type(c_ptr) :: c_array_ptr
    integer(c_int) :: m_size, n_size
    c_array_ptr = rokko_distributed_matrix_get_array_pointer_c(matrix)
    m_size = rokko_distributed_matrix_get_m_size(matrix)
    n_size = rokko_distributed_matrix_get_n_size(matrix)
    call c_f_pointer(c_array_ptr, f_array_ptr, (/m_size,n_size/) )
  end subroutine rokko_distributed_matrix_get_array_pointer
  
end module rokko_distributed_matrix_mod

