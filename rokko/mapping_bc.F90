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

  end interface
  
end module rokko_mapping_bc_mod
