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

module rokko_mapping_bc_mod
  use rokko_grid_mod
  use iso_c_binding
  implicit none

  type, bind(c) :: rokko_mapping_bc
     type(c_ptr) :: ptr
     integer(c_int) :: major
  end type rokko_mapping_bc

  ! generic names
  interface rokko_construct
     procedure rokko_mapping_bc_construct_block_size
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_mapping_bc_destruct
  end interface rokko_destruct

  interface rokko_get_mb
     procedure rokko_mapping_bc_get_mb
  end interface rokko_get_mb

  interface rokko_get_nb
     procedure rokko_mapping_bc_get_nb
  end interface rokko_get_nb

  interface rokko_get_m_global
     procedure rokko_mapping_bc_get_m_global
  end interface rokko_get_m_global

  interface rokko_get_n_global
     procedure rokko_mapping_bc_get_n_global
  end interface rokko_get_n_global

  interface rokko_get_m_local
     procedure rokko_mapping_bc_get_m_local
  end interface rokko_get_m_local

  interface rokko_get_n_local
     procedure rokko_mapping_bc_get_n_local
  end interface rokko_get_n_local

  interface rokko_get_lld
     procedure rokko_mapping_bc_get_lld
  end interface rokko_get_lld

  interface rokko_get_m_size
     procedure rokko_mapping_bc_get_m_size
  end interface rokko_get_m_size

  interface rokko_get_n_size
     procedure rokko_mapping_bc_get_n_size
  end interface rokko_get_n_size

  interface rokko_translate_l2g_row0
     procedure rokko_mapping_bc_translate_l2g_row0
  end interface rokko_translate_l2g_row0

  interface rokko_translate_l2g_col0
     procedure rokko_mapping_bc_translate_l2g_col0
  end interface rokko_translate_l2g_col0

  interface rokko_translate_g2l_row0
     procedure rokko_mapping_bc_translate_g2l_row0
  end interface rokko_translate_g2l_row0

  interface rokko_translate_g2l_col0
     procedure rokko_mapping_bc_translate_g2l_col0
  end interface rokko_translate_g2l_col0

  ! offset by one
  interface rokko_translate_l2g_row
     procedure  rokko_mapping_bc_translate_l2g_row
  end interface rokko_translate_l2g_row

  interface rokko_translate_l2g_col
     procedure  rokko_mapping_bc_translate_l2g_col
  end interface rokko_translate_l2g_col

  interface rokko_translate_g2l_row
     procedure  rokko_mapping_bc_translate_g2l_row
  end interface rokko_translate_g2l_row

  interface rokko_translate_g2l_col
     procedure  rokko_mapping_bc_translate_g2l_col
  end interface rokko_translate_g2l_col

  interface rokko_is_row_major
     procedure rokko_mapping_bc_is_row_major
  end interface rokko_is_row_major

  interface rokko_is_col_major
     procedure rokko_mapping_bc_is_col_major
  end interface rokko_is_col_major

  interface rokko_get_myrank
     procedure rokko_mapping_bc_get_myrank
  end interface rokko_get_myrank

  interface rokko_get_nprocs
     procedure rokko_mapping_bc_get_nprocs
  end interface rokko_get_nprocs

  interface rokko_get_myrow
     procedure rokko_mapping_bc_get_myrow
  end interface rokko_get_myrow

  interface rokko_get_mycol
     procedure rokko_mapping_bc_get_mycol
  end interface rokko_get_mycol

  interface rokko_get_nprow
     procedure rokko_mapping_bc_get_nprow
  end interface rokko_get_nprow

  interface rokko_get_npcol
     procedure rokko_mapping_bc_get_npcol
  end interface rokko_get_npcol

  interface

     subroutine rokko_mapping_bc_construct_block_size(map, global_dim, block_size, grid) &
          bind(c)
       use iso_c_binding
       import rokko_grid, rokko_mapping_bc
       implicit none
       type(rokko_mapping_bc), intent(inout) :: map
       integer(c_int), value, intent(in) :: global_dim
       integer(c_int), value, intent(in) :: block_size
       type(rokko_grid), value, intent(in) :: grid
     end subroutine rokko_mapping_bc_construct_block_size

     subroutine rokko_mapping_bc_destruct(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       type(rokko_mapping_bc), intent(inout) :: map
     end subroutine rokko_mapping_bc_destruct

     function rokko_mapping_bc_get_mb(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_mb
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_mb

     function rokko_mapping_bc_get_nb(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_nb
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_nb

     function rokko_mapping_bc_get_m_local(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_m_local
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_m_local

     function rokko_mapping_bc_get_n_local(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_n_local
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_n_local

     function rokko_mapping_bc_get_m_global(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_m_global
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_m_global

     function rokko_mapping_bc_get_n_global(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_n_global
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_n_global

     function rokko_mapping_bc_get_lld(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_lld
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_lld

     function rokko_mapping_bc_get_m_size(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_m_size
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_m_size

     function rokko_mapping_bc_get_n_size(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_n_size
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_n_size

     function rokko_mapping_bc_translate_l2g_row0(map, local_i) &
          & bind(c,name="rokko_mapping_bc_translate_l2g_row")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_l2g_row0
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int),value,intent(in) :: local_i
     end function rokko_mapping_bc_translate_l2g_row0

     function rokko_mapping_bc_translate_l2g_col0(map, local_j) &
          & bind(c,name="rokko_mapping_bc_translate_l2g_col")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_l2g_col0
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in)::local_j
     end function rokko_mapping_bc_translate_l2g_col0

     function rokko_mapping_bc_translate_g2l_row0(map, global_i) &
          & bind(c,name="rokko_mapping_bc_translate_g2l_row")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_g2l_row0
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in)::global_i
     end function rokko_mapping_bc_translate_g2l_row0

     function rokko_mapping_bc_translate_g2l_col0(map, global_j) &
          & bind(c,name="rokko_mapping_bc_translate_g2l_col")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_g2l_col0
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in) :: global_j
     end function rokko_mapping_bc_translate_g2l_col0

     ! offset by one
     function rokko_mapping_bc_translate_l2g_row(map, local_i) &
          & bind(c,name="rokko_mapping_bc_translate_l2g_row1")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_l2g_row
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int),value,intent(in) :: local_i
     end function rokko_mapping_bc_translate_l2g_row

     function rokko_mapping_bc_translate_l2g_col(map, local_j) &
          & bind(c,name="rokko_mapping_bc_translate_l2g_col1")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_l2g_col
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in)::local_j
     end function rokko_mapping_bc_translate_l2g_col

     function rokko_mapping_bc_translate_g2l_row(map, global_i) &
          & bind(c,name="rokko_mapping_bc_translate_g2l_row1")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_g2l_row
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in)::global_i
     end function rokko_mapping_bc_translate_g2l_row

     function rokko_mapping_bc_translate_g2l_col(map, global_j) &
          & bind(c,name="rokko_mapping_bc_translate_g2l_col1")
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_translate_g2l_col
       type(rokko_mapping_bc), value, intent(in) :: map
       integer(c_int), value, intent(in) :: global_j
     end function rokko_mapping_bc_translate_g2l_col

     function rokko_mapping_bc_is_row_major(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       logical(c_bool) :: rokko_mapping_bc_is_row_major
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_is_row_major

     function rokko_mapping_bc_is_col_major(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       logical(c_bool) :: rokko_mapping_bc_is_col_major
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_is_col_major

     function rokko_mapping_bc_get_myrank(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_myrank
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_myrank

     function rokko_mapping_bc_get_nprocs(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_nprocs
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_nprocs

     function rokko_mapping_bc_get_myrow(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_myrow
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_myrow

     function rokko_mapping_bc_get_mycol(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_mycol
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_mycol

     function rokko_mapping_bc_get_nprow(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_nprow
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_nprow

     function rokko_mapping_bc_get_npcol(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       integer(c_int) :: rokko_mapping_bc_get_npcol
       type(rokko_mapping_bc), value, intent(in) :: map
     end function rokko_mapping_bc_get_npcol

  end interface

end module rokko_mapping_bc_mod
