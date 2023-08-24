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

module rokko_grid_mod
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_grid_col_major = 1, rokko_grid_row_major = 2
  end enum

  type, bind(c) :: rokko_grid
     type(c_ptr) :: ptr
     integer(c_int) :: major
  end type rokko_grid

  ! generic names
  interface rokko_construct
     procedure rokko_grid_construct_f
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_grid_destruct
  end interface rokko_destruct

  interface rokko_get_myrank
     procedure rokko_grid_get_myrank
  end interface rokko_get_myrank

  interface rokko_get_nprocs
     procedure rokko_grid_get_nprocs
  end interface rokko_get_nprocs

  interface rokko_get_myrow
     procedure rokko_grid_get_myrow
  end interface rokko_get_myrow

  interface rokko_get_mycol
     procedure rokko_grid_get_mycol
  end interface rokko_get_mycol

  interface rokko_get_nprow
     procedure rokko_grid_get_nprow
  end interface rokko_get_nprow

  interface rokko_get_npcol
     procedure rokko_grid_get_npcol
  end interface rokko_get_npcol

  interface rokko_is_row_major
     procedure rokko_grid_is_row_major
  end interface rokko_is_row_major

  interface rokko_is_col_major
     procedure rokko_grid_is_col_major
  end interface rokko_is_col_major

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

     function rokko_grid_get_myrow(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       integer(c_int) :: rokko_grid_get_myrow
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_myrow

     function rokko_grid_get_mycol(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       integer(c_int) :: rokko_grid_get_mycol
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_mycol

     function rokko_grid_get_nprow(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       integer(c_int) :: rokko_grid_get_nprow
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_nprow

     function rokko_grid_get_npcol(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       integer(c_int) :: rokko_grid_get_npcol
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_get_npcol

     function rokko_grid_is_row_major(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       logical(c_bool) :: rokko_grid_is_row_major
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_is_row_major

     function rokko_grid_is_col_major(grid) bind(c)
       use iso_c_binding
       import rokko_grid
       implicit none
       logical(c_bool) :: rokko_grid_is_col_major
       type(rokko_grid), value, intent(in) :: grid
     end function rokko_grid_is_col_major

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

end module rokko_grid_mod
