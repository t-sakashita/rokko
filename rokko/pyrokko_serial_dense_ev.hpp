/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/pyrokko_parameters.hpp>
#include <rokko/serial_dense_ev.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/traits/real_t.hpp>

namespace rokko {

class wrap_serial_dense_ev : public serial_dense_ev {
public:
  wrap_serial_dense_ev(std::string const& solver_name) : serial_dense_ev(solver_name) {}

  wrap_serial_dense_ev() = default;

  void initialize() {
    int num = 1;
    char** ptr = NULL;
    serial_dense_ev::initialize(num, ptr);
  }

  template<typename T, int MAJOR>
  wrap_parameters diagonalize(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in,
             Eigen::RefVec<real_t<T>> eigval_in,
             Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> eigvec_in,
			 wrap_parameters const& params) {
    using vector_type = Eigen::Vector<real_t<T>>;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR> eigvec;
    vector_type eigval;

    new (&mat) Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in.data(), mat_in.rows(), mat_in.cols());
    new (&eigval) Eigen::Ref<vector_type>(eigval_in);
    new (&eigvec) Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(eigvec_in.data(), eigvec_in.rows(), eigvec_in.cols());

    const wrap_parameters params_out = serial_dense_ev::diagonalize(mat, eigval, eigvec, parameters(params));

    new (&mat) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
    new (&eigval) vector_type();
    new (&eigvec) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();

    return params_out;
  }

  template<typename T, int MAJOR>
  wrap_parameters diagonalize(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in,
             Eigen::RefVec<real_t<T>> eigval_in,
			 wrap_parameters const& params) {
    using vector_type = Eigen::Vector<real_t<T>>;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat;
    vector_type eigval;

    new (&mat) Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in.data(), mat_in.rows(), mat_in.cols());
    new (&eigval) Eigen::Ref<vector_type>(eigval_in);

    const wrap_parameters params_out = serial_dense_ev::diagonalize(mat, eigval, parameters(params));

    new (&mat) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
    new (&eigval) vector_type();

    return params_out;
  }
};

} // end namespace rokko
