/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ORIGINAL_SOLVER_HPP
#define ROKKO_ORIGINAL_SOLVER_HPP

#include <tuple>

#include <rokko/distributed_vector.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/skel/mapping_1d.hpp>
#include <rokko/parameters.hpp>
#include <rokko/original/lanczos_tridiagonalization.hpp>
#include <rokko/lapack/stebz.hpp>
#include <rokko/lapack/stein.hpp>

namespace rokko {

namespace original {

class solver {
public:
  using value_type = double;

  solver() = default;

  ~solver() = default;

  void initialize(int& argc, char**& argv) {}
  void finalize() {}

  parameters diagonalize(const rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    throw std::invalid_argument("rokko::original::solver::diagonalize() is not implemented.");
  }

  parameters diagonalize(const rokko::distributed_mfree& mat, rokko::parameters const& params) {
    parameters params_out;
    map_ = std::make_shared<rokko::skel::mapping_1d>(mat.get_dim(), mat.get_num_local_rows(), mpi_comm{mat.get_comm()});

    // tridiagonaliation by Lanczos method
    const int dim = mat.get_dim();
    Eigen::VectorXd alpha(dim), beta(dim-1);
    Eigen::MatrixXd u(dim, dim);
    lanczos::tridiagonalization(mat, alpha, beta, u);

    // calculating eigenvalues of tridiagonalized matrix
    eigvals.resize(dim);
    constexpr double abstol = 0.;
    int il = dim, iu = dim;
    int nsplit;
    Eigen::VectorXi iblock(dim), isplit(dim);
    int info = rokko::lapack::stebz('E', il, iu, abstol, alpha, beta, num_conv, nsplit, eigvals, iblock, isplit);

    // calculating eigenvectors of tridiagonalized matrix
    eigvecs.resize(dim, num_conv);
    Eigen::VectorXi ifailv(num_conv);
    info = rokko::lapack::stein(alpha, beta, num_conv, eigvals, eigvecs, iblock, isplit, ifailv);
    // calculating eigenvectors of mat
    eigvecs = u * eigvecs;

    set_output_parameters(params_out);
    return params_out;
  }

  void set_output_parameters(rokko::parameters& params_out) const {
    int num_conv = get_num_conv();
    params_out.set("num_conv", num_conv);
  }

  value_type eigenvalue(int k) const { return eigvals(k); }

  void eigenvector(int k, std::vector<value_type>& vec) const {
    if (vec.size() < map_->get_num_local_rows()) vec.resize(map_->get_num_local_rows());
    eigenvector(k, vec.data());
  }

  void eigenvector(int k, distributed_vector<double>& vec) const {
    vec.initialize(map_->get_dim(), map_->start_row(), map_->end_row());
    eigenvector(k, vec.get_storage());
  }

  void eigenvector(int k, value_type *const vec) const {
    value_type const*const vec_pt = eigvecs.col(k).data();
    std::copy(vec_pt, vec_pt + map_->get_num_local_rows(), vec);
  }

  int get_num_conv() const {
    return num_conv;
  }

private:
  std::shared_ptr<const rokko::skel::mapping_1d> map_;
  Eigen::VectorXd eigvals;
  Eigen::MatrixXd eigvecs;
  int num_conv;
};

} // namespace original

} // namespace rokko

#endif // ROKKO_ORIGINAL_SOLVER_HPP
