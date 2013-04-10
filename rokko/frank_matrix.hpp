#ifndef ROKKO_FRANK_MATRIX_H
#define ROKKO_FRANK_MATRIX_H

#include <rokko/distributed_matrix.hpp>

namespace rokko {

void generate_frank_matrix_local(rokko::distributed_matrix& mat)
{
  for(int local_i=0; local_i<mat.m_local; ++local_i) {
    for(int local_j=0; local_j<mat.n_local; ++local_j) {
      int global_i = mat.translate_l2g_row(local_i);
      int global_j = mat.translate_l2g_col(local_j);
      mat.array[local_i * mat.n_local + local_j] = mat.m_global - max(global_i, global_j);
    }
  }
}

} // namespace rokko

#endif // ROKKO_FRANK_MATRIX_H
