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

module rokko_mapping_bc_mod
!  use rokko_parallel_dense, only : rokko_parallel_dense_ev
!  use rokko_parallel_dense_classes
  use rokko_grid_mod
  use rokko_mapping_bc_type
  use iso_c_binding
  implicit none
  
  interface

     ! subroutine rokko_mapping_bc_construct(map, global_dim, grid, solver) &
     !      bind(c)
     !   use iso_c_binding
     !   import rokko_grid, rokko_mapping_bc, rokko_parallel_dense_ev
     !   implicit none
     !   type(rokko_mapping_bc), intent(inout) :: map
     !   integer(c_int), value, intent(in) :: global_dim
     !   type(rokko_grid), value, intent(in) :: grid
     !   type(rokko_parallel_dense_ev), value, intent(in) :: solver
     ! end subroutine rokko_mapping_bc_construct
     
     subroutine rokko_mapping_bc_destruct(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       type(rokko_mapping_bc), intent(inout) :: map
     end subroutine rokko_mapping_bc_destruct

  end interface
  
end module rokko_mapping_bc_mod
