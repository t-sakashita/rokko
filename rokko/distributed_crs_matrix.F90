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

module rokko_distributed_crs_matrix_mod
  use iso_c_binding
  use rokko_mapping_1d_mod
  implicit none

  type, bind(c) :: rokko_distributed_crs_matrix
     type(c_ptr) :: ptr
  end type rokko_distributed_crs_matrix

  interface

     subroutine rokko_distributed_crs_matrix_construct(matrix, map, num_entries_per_row) &
          & bind(c,name='rokko_distributed_crs_matrix_construct')
       use iso_c_binding
       import rokko_mapping_1d, rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(out) :: matrix
       type(rokko_mapping_1d), value, intent(in) :: map
       integer(c_int), value, intent(in) :: num_entries_per_row
     end subroutine rokko_distributed_crs_matrix_construct

     subroutine rokko_distributed_crs_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_crs_matrix_destruct

     subroutine rokko_distributed_crs_matrix_insert0(matrix, row, col_size, cols, values) &
          & bind(c,name='rokko_distributed_crs_matrix_insert')
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: row, col_size
       integer(c_int), dimension(col_size), intent(in) :: cols
       real(c_double), dimension(col_size), intent(in) :: values
     end subroutine rokko_distributed_crs_matrix_insert0

     subroutine rokko_distributed_crs_matrix_complete(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_crs_matrix_complete

     function rokko_distributed_crs_matrix_start_row0(matrix) &
          & bind(c,name='rokko_distributed_crs_matrix_start_row')
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_start_row0
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end function rokko_distributed_crs_matrix_start_row0

     function rokko_distributed_crs_matrix_end_row0(matrix) &
          & bind(c,name='rokko_distributed_crs_matrix_end_row')
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_end_row0
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end function rokko_distributed_crs_matrix_end_row0

     function rokko_distributed_crs_matrix_num_local_rows(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       integer(c_int) :: rokko_distributed_crs_matrix_num_local_rows
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end function rokko_distributed_crs_matrix_num_local_rows

     subroutine rokko_distributed_crs_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_crs_matrix
       implicit none
       type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_crs_matrix_print
  end interface

  ! generic names
  interface rokko_construct
     procedure rokko_distributed_crs_matrix_construct
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_distributed_crs_matrix_destruct
  end interface rokko_destruct

  interface rokko_insert
     module procedure rokko_distributed_crs_matrix_insert
  end interface rokko_insert

  interface rokko_insert0
     procedure rokko_distributed_crs_matrix_insert0
  end interface rokko_insert0

  interface rokko_complete
     procedure rokko_distributed_crs_matrix_complete
  end interface rokko_complete

  interface rokko_start_row
     module procedure rokko_distributed_crs_matrix_start_row
  end interface rokko_start_row

  interface rokko_end_row
     module procedure rokko_distributed_crs_matrix_end_row
  end interface rokko_end_row

  interface rokko_start_row0
     procedure rokko_distributed_crs_matrix_start_row0
  end interface rokko_start_row0

  interface rokko_end_row0
     procedure rokko_distributed_crs_matrix_end_row0
  end interface rokko_end_row0

  interface rokko_num_local_rows
     procedure rokko_distributed_crs_matrix_num_local_rows
  end interface rokko_num_local_rows

  interface rokko_print
     procedure rokko_distributed_crs_matrix_print
  end interface rokko_print

contains

  function rokko_distributed_crs_matrix_start_row(matrix) result(ind)
    integer :: ind
    type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
    ind = rokko_distributed_crs_matrix_start_row0(matrix) + 1
  end function rokko_distributed_crs_matrix_start_row

  function rokko_distributed_crs_matrix_end_row(matrix) result(ind)
    integer :: ind
    type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
    ind = rokko_distributed_crs_matrix_end_row0(matrix)
  end function rokko_distributed_crs_matrix_end_row

  subroutine rokko_distributed_crs_matrix_insert(matrix, row, col_size, cols, values)
    type(rokko_distributed_crs_matrix), value, intent(in) :: matrix
    integer, intent(in) :: row, col_size
    integer, dimension(col_size), intent(in) :: cols
    integer, dimension(col_size) :: cols2
    double precision, dimension(col_size), intent(in) :: values
    integer :: i

    do i=1, col_size
       cols2(i) = cols(i) - 1
    enddo
    call rokko_distributed_crs_matrix_insert0(matrix, row-1, col_size, cols2, values)
  end subroutine rokko_distributed_crs_matrix_insert

end module rokko_distributed_crs_matrix_mod
