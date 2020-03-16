!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module rokko_mapping_1d_mod
  use mpi
  use iso_c_binding
  implicit none

  type, bind(c) :: rokko_mapping_1d
     type(c_ptr) :: ptr
  end type rokko_mapping_1d

  ! generic names
  interface rokko_construct
     procedure rokko_mapping_1d_construct
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_mapping_1d_destruct
  end interface rokko_destruct

  interface rokko_get_dim
     procedure rokko_mapping_1d_get_dim
  end interface rokko_get_dim

  interface rokko_num_local_rows
     procedure rokko_mapping_1d_num_local_rows
  end interface rokko_num_local_rows

  interface rokko_start_row
     module procedure rokko_mapping_1d_start_row
  end interface rokko_start_row

  interface rokko_end_row
     module procedure rokko_mapping_1d_end_row
  end interface rokko_end_row

  interface rokko_start_row0
     procedure rokko_mapping_1d_start_row0
  end interface rokko_start_row0

  interface rokko_end_row0
     procedure rokko_mapping_1d_end_row0
  end interface rokko_end_row0

  interface rokko_get_comm
     procedure rokko_mapping_1d_get_comm
  end interface rokko_get_comm

  interface

     subroutine rokko_mapping_1d_construct(map, dim, comm) &
          & bind(c,name='rokko_mapping_1d_construct_f')
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       type(rokko_mapping_1d), intent(inout) :: map
       integer(c_int), value, intent(in) :: dim
       integer, value, intent(in) :: comm
     end subroutine rokko_mapping_1d_construct

     subroutine rokko_mapping_1d_destruct(map) bind(c)
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       type(rokko_mapping_1d), intent(inout) :: map
     end subroutine rokko_mapping_1d_destruct

     function rokko_mapping_1d_get_dim(map) bind(c)
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       integer(c_int) :: rokko_mapping_1d_get_dim
       type(rokko_mapping_1d), value, intent(in) :: map
     end function rokko_mapping_1d_get_dim

     function rokko_mapping_1d_num_local_rows(map) bind(c)
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       integer(c_int) :: rokko_mapping_1d_num_local_rows
       type(rokko_mapping_1d), value, intent(in) :: map
     end function rokko_mapping_1d_num_local_rows

     function rokko_mapping_1d_start_row0(map) &
          & bind(c,name='rokko_mapping_1d_start_row')
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       integer(c_int) :: rokko_mapping_1d_start_row0
       type(rokko_mapping_1d), value, intent(in) :: map
     end function rokko_mapping_1d_start_row0

     function rokko_mapping_1d_end_row0(map) &
          & bind(c,name='rokko_mapping_1d_end_row')
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       integer(c_int) :: rokko_mapping_1d_end_row0
       type(rokko_mapping_1d), value, intent(in) :: map
     end function rokko_mapping_1d_end_row0

     function rokko_mapping_1d_get_comm(map) &
          & bind(c,name='rokko_mapping_1d_get_comm_f')
       use iso_c_binding
       import rokko_mapping_1d
       implicit none
       integer(c_int) :: rokko_mapping_1d_get_comm
       type(rokko_mapping_1d), value, intent(in) :: map
     end function rokko_mapping_1d_get_comm

  end interface

contains

  function rokko_mapping_1d_start_row(matrix) result(ind)
    integer :: ind
    type(rokko_mapping_1d), value, intent(in) :: matrix
    ind = rokko_mapping_1d_start_row0(matrix) + 1
  end function rokko_mapping_1d_start_row

  function rokko_mapping_1d_end_row(matrix) result(ind)
    integer :: ind
    type(rokko_mapping_1d), value, intent(in) :: matrix
    ind = rokko_mapping_1d_end_row0(matrix)
  end function rokko_mapping_1d_end_row

end module rokko_mapping_1d_mod
