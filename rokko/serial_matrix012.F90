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

module rokko_serial_matrix012_mod
  use rokko_eigen_matrix_mod
  implicit none

  interface rokko_matrix012_generate
     procedure rokko_matrix012_generate_eigen_matrix
  end interface rokko_matrix012_generate

  interface

     subroutine rokko_matrix012_generate_eigen_matrix(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
     end subroutine rokko_matrix012_generate_eigen_matrix

  end interface

contains

end module rokko_serial_matrix012_mod
