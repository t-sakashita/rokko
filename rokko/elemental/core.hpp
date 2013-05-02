#ifndef ROKKO_CORE_ELEMENTAL_H
#define ROKKO_CORE_ELEMENTAL_H

#include <rokko/elemental/elemental.hpp>

//#include <boost/noncopyable.hpp>

#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>

namespace rokko {


class solver_elemental
{
public:
  void initialize(int& argc, char**& argv)
  {
    elem::Initialize(argc, argv);
  }

  void finalize()
  {
    elem::Finalize();
  }

  void optimized_grid_size()
  {
  }

  void optimized_matrix_size()
  {
  }

  template<typename MATRIX_MAJOR>
  void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs);

};


} // namespace rokko

#endif // ROKKO_CORE_ELEMENTAL_H
