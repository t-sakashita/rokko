!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine rokko_distributed_matrix_generate_array(matrix, array)
  use rokko, only: rokko_distributed_matrix, rokko_distributed_matrix_get_m_local, &
       rokko_distributed_matrix_get_n_local, rokko_distributed_matrix_translate_l2g_row, &
       rokko_distributed_matrix_translate_l2g_col, rokko_distributed_matrix_set_local
  implicit none
  type(rokko_distributed_matrix), intent(out) :: matrix
  real(8), intent(in) :: array(:,:)
  integer :: m_local, n_local, local_i, local_j, global_i, global_j
  m_local = rokko_distributed_matrix_get_m_local(matrix)
  n_local = rokko_distributed_matrix_get_n_local(matrix)
  do local_i = 0, m_local - 1 
     do local_j = 0, n_local - 1
        global_i = rokko_distributed_matrix_translate_l2g_row(matrix, local_i)
        global_j = rokko_distributed_matrix_translate_l2g_col(matrix, local_j)
        call rokko_distributed_matrix_set_local(matrix, local_i, local_j, &
             array(global_i + 1, global_j + 1))
     enddo
  enddo
end subroutine rokko_distributed_matrix_generate_array
