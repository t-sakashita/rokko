!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>,
!                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module rokko_frank_matrix
  use iso_c_binding
  implicit none
  
  interface

     subroutine rokko_frank_matrix_generate_localized_matrix(matrix) bind(c)
       use rokko
       implicit none
       type(rokko_localized_matrix), intent(inout) :: matrix
     end subroutine rokko_frank_matrix_generate_localized_matrix

     subroutine rokko_frank_matrix_generate_distributed_matrix(matrix) bind(c)
       use rokko
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_frank_matrix_generate_distributed_matrix

  end interface
end module rokko_frank_matrix
