#ifndef ROKKO_EIGEN_SX_CORE_HPP
#define ROKKO_EIGEN_SX_CORE_HPP

#include <rokko/eigen_sx/eigen_sx.hpp>
#include <rokko/eigen_sx/diagonalize.hpp>

namespace rokko {
namespace eigen_sx {

class solver {
public:
  void initialize(int& argc, char**& argv) {}

  void finalize() {}

  template<typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    int n = mat.get_m_global();
    int nx = (n-1) / mat.get_nprow() + 1;
    int i1 = 6, i2 = 16 * 2, i3 = 16 * 4, nm;
    CSTAB_get_optdim(nx, i1, i2, i3, nm);  // return an optimized (possiblly) leading dimension of local block-cyclic matrix to nm.
    //int para_int = 0;   eigen_free_wrapper(para_int);

    int NB  = 64 + 32;
    int nmz = (n-1) / mat.get_nprow() + 1;
    nmz = ((nmz - 1) / NB + 1) * NB + 1;
    int nmw = (n-1) / mat.get_npcol() + 1;
    nmw = ((nmw - 1) / NB + 1) * NB + 1;
    ///int nme = ((n - 1) / 2 + 1) * 2;
    int nh = (n - 1) / 4 + 1;
    i1 = 4;
    int nnh;
    CSTAB_get_optdim(nh, i1, i2, i3, nnh);
    nnh = 4 * nnh;
    int n1x = (n-1) / mat.get_nprocs() + 1;

#ifndef NDEBUG
    std::cout << "nm=" << nm << std::endl;
    std::cout << "nmz=" << nmz << std::endl;
    std::cout << "nmw=" << nmw << std::endl;
    //std::cout << "nme=" << nme << std::endl;
    std::cout << "nnh=" << nnh << std::endl;
    std::cout << "n1x=" << n1x << std::endl;
    //std::cout << "length_array=" << mat.get_length_array() << std::endl;
#endif

    // calculate sizes of my proc's local part of distributed matrix
    mat.set_lld(nm);
    mat.set_length_array(std::max(std::max(std::max(nmz, nm), nnh) * nmw, n*n1x));
    mat.set_block_size(1, 1);
    int m_local = mat.calculate_row_size();
    int n_local = mat.calculate_col_size();
    mat.set_local_size(m_local, n_local);
  }

  template<typename MATRIX_MAJOR>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
    rokko::eigen_sx::diagonalize(mat, eigvals, eigvecs, timer_in);
  }
};

} // namespace eigen_sx
} // namespace rokko

#endif // ROKKO_EIGEN_SX_CORE_HPP

