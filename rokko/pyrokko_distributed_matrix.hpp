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

#ifndef PYROKKO_DISTRIBUTED_MATRIX_HPP
#define PYROKKO_DISTRIBUTED_MATRIX_HPP

#include <pybind11/pybind11.h>

#include <rokko/pyrokko_mapping_bc.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/utility/tuple_to_array.hpp>

#include <memory>

namespace rokko {

class base_distributed_matrix {
public:
  virtual ~base_distributed_matrix() = default;
};

template <typename T, typename MATRIX_MAJOR>
class wrap_distributed_matrix : public distributed_matrix<T,MATRIX_MAJOR>, public base_distributed_matrix {
public:
  using distributed_matrix<T,MATRIX_MAJOR>::get_global_size;
  using distributed_matrix<T,MATRIX_MAJOR>::get_local_size;
  using distributed_matrix<T,MATRIX_MAJOR>::get_block_size;
  using distributed_matrix<T,MATRIX_MAJOR>::get_m_local;
  using distributed_matrix<T,MATRIX_MAJOR>::get_n_local;
  using distributed_matrix<T,MATRIX_MAJOR>::get_lld;
  using distributed_matrix<T,MATRIX_MAJOR>::get_array_pointer;
  using distributed_matrix<T,MATRIX_MAJOR>::is_col_major;

  wrap_distributed_matrix(wrap_mapping_bc<MATRIX_MAJOR> const& map) : distributed_matrix<T,MATRIX_MAJOR>(map) {}

  wrap_distributed_matrix() = default;
  
  std::tuple<int,int> get_block_shape() const {
    return get_block_size();
  }
  
  bool has_global_indices(std::tuple<int,int> const& global) const {
    return has_global_indices(to_array(global));
  }
  
  std::tuple<int,int> get_global_shape() const {
    return get_global_size();
  }
  
  std::tuple<int,int> get_local_shape() const {
    return get_local_size();
  }

  std::tuple<int,int> translate_l2g(std::tuple<int,int> const& local) const {
    return translate_l2g(to_array(local));
  }

  std::tuple<int,int> translate_g2l(std::tuple<int,int> const& global) const {
    return translate_g2l(to_array(global));
  }

  void print() const {
    distributed_matrix<T,MATRIX_MAJOR>::print();
  }

  std::string get_major_string() const {
    return is_col_major() ? "col" : "row";
  }

  wrap_mapping_bc<MATRIX_MAJOR> get_map() const {
    return distributed_matrix<T,MATRIX_MAJOR>::get_mapping();
  }

  py::array_t<T> get_ndarray() {
    if (is_col_major())
      return py::array_t<T>(get_local_size(), {sizeof(T), sizeof(T)*get_lld()}, get_array_pointer(), py::cast(*this));
    else
     return py::array_t<T>(get_local_size(), {sizeof(T)*get_lld(), sizeof(T)}, get_array_pointer(), py::cast(*this));
  }

  auto get_eigen_map() {
    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, rokko::eigen3_major<MATRIX_MAJOR>>,0,Eigen::OuterStride<>>(get_array_pointer(), get_m_local(), get_n_local(), Eigen::OuterStride<Eigen::Dynamic>(get_lld()));
  }

  void set_ndarray(py::array_t<T> const& mat) {
    py::array_t<T> array;
    if (is_col_major()) {
      array = py::array_t<T>(get_local_size(), {sizeof(T), sizeof(T)*get_lld()}, get_array_pointer(), py::cast(*this));
    } else {
      array = py::array_t<T>(get_local_size(), {sizeof(T)*get_lld(), sizeof(T)}, get_array_pointer(), py::cast(*this));
    }

    auto r = array.template mutable_unchecked<2>();
    for (auto i = 0; i < r.shape(0); ++i) {
      for (auto j = 0; j < r.shape(1); ++j) {
        r(i, j) = *mat.data(i, j);
      }
    }
  }

  void set_col_major_matrix(Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> mat) {
    if (!is_col_major())
      throw std::invalid_argument("Cannot set col-major ndarray to row-major distributed_matrix");

    get_eigen_map() = static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(mat);
  }

  void set_row_major_matrix(Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mat) {
    if (is_col_major())
      throw std::invalid_argument("Cannot set row-major ndarray to col-major distributed_matrix");

    get_eigen_map() = static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(mat);
  }
};


template<typename MATRIX_MAJOR>
std::shared_ptr<base_distributed_matrix> create_distributed_matrix(wrap_mapping_bc<MATRIX_MAJOR> const& map, std::string const& dtype = "double") {
  if (dtype == "double")
    return std::make_shared<wrap_distributed_matrix<double,MATRIX_MAJOR>>(map);
  else if (dtype == "float")
    return std::make_shared<wrap_distributed_matrix<float,MATRIX_MAJOR>>(map);
  else if (dtype == "cdouble")
    return std::make_shared<wrap_distributed_matrix<std::complex<double>,MATRIX_MAJOR>>(map);
  else if (dtype == "cfloat")
    return std::make_shared<wrap_distributed_matrix<std::complex<float>,MATRIX_MAJOR>>(map);
  else
    throw std::invalid_argument("create_distributed_matrix : unknwon type '" + dtype + "'");
}

template<typename T, typename MATRIX_MAJOR>
void pyrokko_product(T alpha,
                     wrap_distributed_matrix<T, matrix_col_major> const& matA, bool transA,
                     wrap_distributed_matrix<T, matrix_col_major> const& matB, bool transB,
                     T beta,
                     wrap_distributed_matrix<T, matrix_col_major>& matC) {
  product(alpha, matA, transA, matB, transB, beta, matC);
}

} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_MATRIX_HPP
