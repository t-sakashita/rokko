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

#pragma once

#include <tuple>

#include <rokko/distributed_vector.hpp>
#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/skel/mapping_1d.hpp>
#include <rokko/parameters.hpp>
#include <rokko/original/distributed_crs_matrix.hpp>
#include <rokko/original/lanczos_tridiagonalization.hpp>
#include <rokko/lapack/diagonalize_tridiagonal_matrix_bisection.hpp>

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
    if (mat.get_solver_name() != "original") {
      throw std::invalid_argument("rokko::original::solver::diagonalize() : " + mat.get_solver_name() + "'s distributed_crs_matrix is given.");
    }
    return diagonalize(*std::static_pointer_cast<const rokko::original::distributed_crs_matrix>(mat.get_ptr()), params);
  }

  parameters diagonalize(const rokko::original::distributed_crs_matrix& mat, rokko::parameters const& params) {
    map_ = std::make_shared<rokko::skel::mapping_1d>(mat.get_map());

    // tridiagonaliation by Lanczos method
    const auto dim = mat.get_dim();
    Eigen::VectorXd alpha(dim), beta(dim-1);
    Eigen::MatrixXd u(dim, dim);
    lanczos::tridiagonalization(*mat.get_matrix(), alpha, beta, u);

    // calculating eigenvalues & eigenvectors of tridiagonalized matrix
    auto params_out = rokko::lapack::diagonalize_bisection(alpha, beta, eigvals, eigvecs, params);
    num_conv = params_out.get<int>("num_conv");

    // calculating eigenvectors of mat
    eigvecs = u * eigvecs;

    return params_out;
  }

  parameters diagonalize(const rokko::distributed_mfree& mat, rokko::parameters const& params) {
    map_ = std::make_shared<rokko::skel::mapping_1d>(mat.get_dim(), mat.get_num_local_rows(), mpi_comm{mat.get_comm()});

    // tridiagonaliation by Lanczos method
    const auto dim = mat.get_dim();
    Eigen::VectorXd alpha(dim), beta(dim-1);
    Eigen::MatrixXd u(dim, dim);
    lanczos::tridiagonalization(mat, alpha, beta, u);

    // calculating eigenvalues & eigenvectors of tridiagonalized matrix
    auto params_out = rokko::lapack::diagonalize_bisection(alpha, beta, eigvals, eigvecs, params);
    num_conv = params_out.get<int>("num_conv");

    // calculating eigenvectors of mat
    eigvecs = u * eigvecs;

    return params_out;
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
