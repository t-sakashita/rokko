#ifndef ROKKO_ELEMENTAL_CORE_H
#define ROKKO_ELEMENTAL_CORE_H

#include <elemental.hpp>
#include <rokko/elemental/diagonalize.hpp>

namespace rokko {
namespace elemental {

class solver {
public:
  void initialize(int& argc, char**& argv) { elem::Initialize(argc, argv); }

  void finalize() { elem::Finalize(); }

  void optimized_grid_size() {}

  void optimized_matrix_size() {}

  template<typename MATRIX_MAJOR>
  void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals,
    rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    rokko::elemental::diagonalize(mat, eigvals, eigvecs);
  }
};

} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_CORE_H
