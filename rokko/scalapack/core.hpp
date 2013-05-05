#ifndef ROKKO_SCALAPACK_CORE_HPP
#define ROKKO_SCALAPACK_CORE_HPP

#include <rokko/scalapack/scalapack.hpp>
#include <rokko/scalapack/diagonalize.hpp>

namespace rokko {
namespace scalapack {

class solver {
public:
  void initialize(int& argc, char**& argv) { }

  void finalize() { }

  void optimized_grid_size() {}

  void optimized_matrix_size(int dim, int nprow, int npcol, int& mb, int& nb, int& lld, int& len_array) {}

  template<typename MATRIX_MAJOR>
  void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals,
                   rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    rokko::scalapack::diagonalize(mat, eigvals, eigvecs);
  }
};

} // namespace sclapack
} // namespace rokko


#endif // ROKKO_SCALAPACK_CORE_HPP

