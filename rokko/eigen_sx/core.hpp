#ifndef ROKKO_EIGEN_SX_CORE_HPP
#define ROKKO_EIGEN_SX_CORE_HPP

#include <rokko/eigen_sx/eigen_sx.hpp>

#include <boost/noncopyable.hpp>

#include <rokko/core.hpp>

namespace rokko {

template<>
void Initialize<rokko::eigen_sx>(int& argc, char**& argv)
{
  int size_of_col_local, size_of_row_local;
  int ndims = 2;
  eigen_init_wrapper(ndims, size_of_col_local, size_of_row_local);

  int nprow = cycl2d_.size_of_row;
  int npcol = cycl2d_.size_of_col;
  int myrow = cycl2d_.my_row;
  int mycol = cycl2d_.my_col;
  //cout << "NPROW=" << cycl2d_.size_of_row << "  NPCOL=" << cycl2d_.size_of_col  << endl;
}

template<>
void Finalize<rokko::eigen_sx>()
{
}

} // namespace rokko

#endif // ROKKO_EIGEN_SX_CORE_HPP

