#ifndef ROKKO_CHECK_ORTHOGONALITY_HPP
#define ROKKO_CHECK_ORTHOGONALITY_HPP

#include <rokko/distributed_matrix.hpp>

namespace rokko {

template<typename MATRIX_MAJOR>
int check_orthogonality(rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs)
{

  rokko::distributed_matrix mat2(mat.get_m_global(), mat.get_n_global(), mat.g);
  product(mat, false, mat, true, 1, 0, mat2);
  std::cout << "check_orthogonality";
  mat2.print();

  return 0;
}

} // namespace rokko

#endif // ROKKO_CHECK_ORTHOGONALITY_HPP
