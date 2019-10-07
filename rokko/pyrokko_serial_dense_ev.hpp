/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_SERIAL_DENSE_EV_HPP
#define PYROKKO_SERIAL_DENSE_EV_HPP

#include <rokko/pyrokko_localized_matrix.hpp>
#include <rokko/pyrokko_parameters.hpp>
#include <rokko/serial_dense_ev.hpp>

namespace rokko {

class wrap_serial_dense_ev : public serial_dense_ev {
public:
  wrap_serial_dense_ev(std::string const& solver_name) : serial_dense_ev(solver_name) {}
  
  wrap_serial_dense_ev() {}

  void initialize() {
    int num = 1;
    char** ptr = NULL;
    serial_dense_ev::initialize(num, ptr);
  }

  template<int MAJOR>
  wrap_parameters diagonalize(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in,
             Eigen::Ref<Eigen::VectorXd> eigval_in,
             Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> eigvec_in,
			 wrap_parameters const& params) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> eigvec;
    Eigen::VectorXd eigval;

    new (&mat) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in.data(), mat_in.rows(), mat_in.cols());
    new (&eigval) Eigen::Ref<Eigen::VectorXd>(eigval_in);
    new (&eigvec) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(eigvec_in.data(), eigvec_in.rows(), eigvec_in.cols());

    wrap_parameters params_out = serial_dense_ev::diagonalize(mat, eigval, eigvec, parameters(params));

    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
    new (&eigval) Eigen::VectorXd();
    new (&eigvec) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();

    return params_out;
  }

  template<int MAJOR>
  wrap_parameters diagonalize(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in,
             Eigen::Ref<Eigen::Vector<double, Eigen::Dynamic>> eigval_in,
			 wrap_parameters const& params) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat;
    Eigen::Vector<double,Eigen::Dynamic> eigval;

    new (&mat) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in.data(), mat_in.rows(), mat_in.cols());
    new (&eigval) Eigen::Ref<Eigen::VectorXd>(eigval_in);

    wrap_parameters params_out = serial_dense_ev::diagonalize(mat, eigval, parameters(params));

    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
    new (&eigval) Eigen::Vector<double, Eigen::Dynamic>();

    return params_out;
  }
};

} // end namespace rokko

#endif // PYROKKO_SERIAL_DENSE_EV_HPP
