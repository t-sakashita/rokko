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

module rokko_localized_matrix_mod
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
  end enum

  type, bind(c) :: rokko_localized_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_localized_matrix

  ! generic names
  interface construct
     procedure rokko_localized_matrix_construct
  end interface construct

  interface destruct
     procedure rokko_localized_matrix_destruct
  end interface destruct

  interface print
     procedure rokko_localized_matrix_print
  end interface print
  
  interface
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

     subroutine rokko_localized_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), value, intent(in) :: matrix
     end subroutine rokko_localized_matrix_print
  end interface

  !
  ! rokko_frank_matrix for localized_matrix
  !

  interface frank_matrix_generate
     procedure rokko_frank_matrix_generate_localized_matrix
  end interface frank_matrix_generate
  
  interface
     subroutine rokko_frank_matrix_generate_localized_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_localized_matrix
       implicit none
       type(rokko_localized_matrix), value, intent(in) :: matrix
     end subroutine rokko_frank_matrix_generate_localized_matrix
  end interface

end module rokko_localized_matrix_mod
