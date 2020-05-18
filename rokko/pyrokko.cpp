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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <rokko/utility/tuple_to_array.hpp>

#include <rokko/eigen3.hpp>
#include <rokko/pyrokko_serial_dense_ev.hpp>
#include <rokko/utility/pyrokko_generate_matrix.hpp>

#include <rokko/pyrokko_parameters.hpp>
#include <rokko/pyrokko_grid.hpp>
#include <rokko/pyrokko_mapping_bc.hpp>
#include <rokko/pyrokko_distributed_matrix.hpp>
#include <rokko/pyrokko_parallel_dense_ev.hpp>
#include <rokko/pyrokko_collective.hpp>

#include <rokko/utility/pyrokko_minij_matrix.hpp>
#include <rokko/utility/pyrokko_frank_matrix.hpp>
#include <rokko/utility/pyrokko_helmert_matrix.hpp>
#include <rokko/utility/pyrokko_laplacian_matrix.hpp>
#include <rokko/utility/pyrokko_matrix012.hpp>

#include <rokko/pyrokko_parallel_sparse_ev.hpp>
#include <rokko/pyrokko_mapping_1d.hpp>
#include <rokko/pyrokko_distributed_crs_matrix.hpp>
#include <rokko/pyrokko_distributed_mfree.hpp>

#include <rokko/utility/pyrokko_sort_eigenpairs.hpp>
#include <rokko/utility/pyrokko_solver_name.hpp>


namespace rokko {

namespace py = pybind11;

PYBIND11_MODULE(pyrokko, m) {
  py::class_<wrap_parameters>(m, "parameters")
    .def(py::init<>())
    .def("clear", py::overload_cast<>(&parameters::clear))
    .def("clear", py::overload_cast<parameters::key_type const&>(&parameters::clear))
    .def("set", static_cast<void (parameters::*)(parameters::key_type const&, const char *)>(&parameters::set))
    .def("set", &parameters::set<double>, py::arg(), py::arg("value").noconvert())
    .def("set", &parameters::set<bool>, py::arg(), py::arg("value").noconvert()) // needed to place before <int>
    .def("set", &parameters::set<int>, py::arg(), py::arg("value").noconvert())
    //.def("set", py::overload_cast<parameters::key_type const&, double const&>(&parameters::set<double>))
    //.def("set", py::overload_cast<parameters::key_type const&, bool const&>(&parameters::set<bool>))
    //.def("set", py::overload_cast<parameters::key_type const&, int const&>(&parameters::set<int>))
    .def("defined", &parameters::defined)
    .def_property_readonly("keys", &parameters::keys)
    .def("get", &wrap_parameters::python_get)
    //.def("get", &parameters::get<int>)
    //.def("get", &parameters::get<double>)
    //.def("get", &parameters::get<bool>)
    .def("get_string", &parameters::get_string)
    .def("get_bool", &parameters::get_bool)
    .def_property_readonly("dict", &wrap_parameters::dict);

  py::enum_<matrix_major_enum>(m, "matrix_major")
    .value("row", matrix_major_enum::row)
    .value("col", matrix_major_enum::col);

  py::class_<wrap_serial_dense_ev>(m, "serial_dense_ev")
    .def(py::init<std::string const&>())
    .def(py::init<>())
    //.def("initialize", py::overload_cast<int&, char**&>(&serial_dense_ev::initialize))
    .def("finalize", &serial_dense_ev::finalize)
    // for double
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<double,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<double,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<double,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<double,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    // for complex double
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<double>,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<double>,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<double>,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<double>,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    // for float
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXf>,Eigen::Ref<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<float,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXf>,Eigen::Ref<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<float,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXf>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<float,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXf>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<float,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    // for complex float
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXf>,Eigen::Ref<Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<float>,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXf>,Eigen::Ref<Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<float>,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXf>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<float>,Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXf>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<std::complex<float>,Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    .def_property_readonly_static("solvers", &serial_dense_ev::solvers)
    .def_property_readonly_static("default_solver", &serial_dense_ev::default_solver);


  // For grid
  py::class_<grid_row_major_t>(m,"grid_row_major")
    .def(py::init<>());
  m.attr("grid_row_major") = py::cast(grid_row_major);
  
  py::class_<grid_col_major_t>(m,"grid_col_major")
    .def(py::init<>());
  m.attr("grid_col_major") = py::cast(grid_col_major);
  
  py::class_<wrap_grid>(m, "grid")
    .def(py::init<>())
    .def(py::init<grid_row_major_t>())
    .def(py::init<grid_col_major_t>())
    .def(py::init<pybind11::handle const&, std::tuple<int,int> const&, grid_row_major_t>(), py::arg("comm"), py::arg("size"), py::arg("major") = grid_row_major)
    .def(py::init<pybind11::handle const&, std::tuple<int,int> const&, grid_col_major_t>(), py::arg("comm"), py::arg("size"), py::arg("major") = grid_col_major)
    .def(py::init<pybind11::handle const&, grid_row_major_t>(), py::arg("comm"), py::arg("major") = grid_row_major)
    .def(py::init<pybind11::handle const&, grid_col_major_t>(), py::arg("comm"), py::arg("major") = grid_col_major)
    .def("calculate_rank_form_coords", &wrap_grid::calculate_rank_form_coords)
    .def("get_comm", &wrap_grid::get_comm)
    .def_property_readonly("comm", &wrap_grid::get_comm)
    .def_property_readonly("nprocs", &wrap_grid::get_nprocs)
    .def_property_readonly("nprow", &wrap_grid::get_nprow)
    .def_property_readonly("npcol", &wrap_grid::get_npcol)
    .def_property_readonly("myrank", &wrap_grid::get_myrank)
    .def_property_readonly("myrow", &wrap_grid::get_myrow)
    .def_property_readonly("mycol", &wrap_grid::get_mycol)
    .def_property_readonly("shape", &wrap_grid::get_shape)
    .def_property_readonly("mine", &wrap_grid::get_mine)
    .def_property_readonly("major", &wrap_grid::get_major_string);

  py::class_<base_mapping_bc, std::shared_ptr<base_mapping_bc>>(m, "base_mapping_bc");

  py::class_<wrap_mapping_bc<matrix_col_major>, base_mapping_bc, std::shared_ptr<wrap_mapping_bc<matrix_col_major>>>(m, "mapping_bc_col")
    .def(py::init<int, int, wrap_grid>(), py::arg("global_dim"), py::arg("block_size"), py::arg("grid"))
    .def(py::init<int, int, int, wrap_grid>(), py::arg("global_dim"), py::arg("block_size"), py::arg("lld"), py::arg("grid"))
    .def(py::init<std::tuple<int,int> const&, std::tuple<int,int> const&, wrap_grid>(), py::arg("global_size"), py::arg("block_size"), py::arg("grid"))
    .def("get_mb", &wrap_mapping_bc<matrix_col_major>::get_mb)
    .def("get_nb", &wrap_mapping_bc<matrix_col_major>::get_nb)
    .def_property_readonly("block_shape", &wrap_mapping_bc<matrix_col_major>::get_block_shape)
    .def_property_readonly("m_global", &wrap_mapping_bc<matrix_col_major>::get_m_global)
    .def_property_readonly("n_global", &wrap_mapping_bc<matrix_col_major>::get_n_global)
    .def_property_readonly("global_shape", &wrap_mapping_bc<matrix_col_major>::get_global_shape)
    .def_property_readonly("m_local", &wrap_mapping_bc<matrix_col_major>::get_m_local)
    .def_property_readonly("n_local", &wrap_mapping_bc<matrix_col_major>::get_n_local)
    .def_property_readonly("local_shape", &wrap_mapping_bc<matrix_col_major>::get_local_shape)
    .def("translate_l2g_row", &wrap_mapping_bc<matrix_col_major>::translate_l2g_row)
    .def("translate_l2g_col", &wrap_mapping_bc<matrix_col_major>::translate_l2g_col)
    .def("translate_l2g", &wrap_mapping_bc<matrix_col_major>::translate_l2g) // tuple
    .def("translate_g2l_row", &wrap_mapping_bc<matrix_col_major>::translate_g2l_row)
    .def("translate_g2l_col", &wrap_mapping_bc<matrix_col_major>::translate_g2l_col)
    .def("translate_g2l", &wrap_mapping_bc<matrix_col_major>::translate_g2l) // tuple
    .def("has_global_row_index", &wrap_mapping_bc<matrix_col_major>::has_global_row_index)
    .def("has_global_col_index", &wrap_mapping_bc<matrix_col_major>::has_global_col_index)
    .def("has_global_indices", &wrap_mapping_bc<matrix_col_major>::has_global_indices)
    .def("get_grid", &wrap_mapping_bc<matrix_col_major>::get_grid)
    .def_property_readonly("major", &wrap_mapping_bc<matrix_col_major>::get_major_string);

  m.def("mapping_bc", py::overload_cast<std::tuple<int,int> const&, std::tuple<int,int> const&, wrap_grid const&, matrix_major_enum const&>(&create_mapping_bc));
  m.def("mapping_bc", py::overload_cast<int, int, int, wrap_grid const&, matrix_major_enum const&>(&create_mapping_bc));
  m.def("mapping_bc", py::overload_cast<int, int, wrap_grid const&, matrix_major_enum const&>(&create_mapping_bc));

  py::class_<wrap_distributed_matrix<matrix_col_major>>(m, "distributed_matrix")
    .def(py::init<wrap_mapping_bc<matrix_col_major>>())
    .def_property_readonly("mb", &wrap_distributed_matrix<matrix_col_major>::get_mb)
    .def_property_readonly("nb", &wrap_distributed_matrix<matrix_col_major>::get_nb)
    .def_property_readonly("block_shape", &wrap_distributed_matrix<matrix_col_major>::get_block_shape)
    .def_property_readonly("m_global", &wrap_distributed_matrix<matrix_col_major>::get_m_global)
    .def_property_readonly("n_global", &wrap_distributed_matrix<matrix_col_major>::get_n_global)
    .def_property_readonly("m_local", &wrap_distributed_matrix<matrix_col_major>::get_m_local)
    .def_property_readonly("n_local", &wrap_distributed_matrix<matrix_col_major>::get_n_local)
    .def_property_readonly("global_shape", &wrap_distributed_matrix<matrix_col_major>::get_global_shape)
    .def_property_readonly("local_shape", &wrap_distributed_matrix<matrix_col_major>::get_local_shape)
    .def("translate_l2g_row", &wrap_distributed_matrix<matrix_col_major>::translate_l2g_row)
    .def("translate_l2g_col", &wrap_distributed_matrix<matrix_col_major>::translate_l2g_col)
    .def("translate_l2g", &wrap_distributed_matrix<matrix_col_major>::translate_l2g) // tuple
    .def("translate_g2l_row", &wrap_distributed_matrix<matrix_col_major>::translate_g2l_row)
    .def("translate_g2l_col", &wrap_distributed_matrix<matrix_col_major>::translate_g2l_col)
    .def("translate_g2l", &wrap_distributed_matrix<matrix_col_major>::translate_g2l) // tuple
    .def("has_global_row_index", &wrap_distributed_matrix<matrix_col_major>::has_global_row_index)
    .def("has_global_col_index", &wrap_distributed_matrix<matrix_col_major>::has_global_col_index)
    .def("has_global_indices", &wrap_distributed_matrix<matrix_col_major>::has_global_indices)
    .def("get_local", &wrap_distributed_matrix<matrix_col_major>::get_local)
    .def("get_global", &wrap_distributed_matrix<matrix_col_major>::get_global)
    .def("set_local", &wrap_distributed_matrix<matrix_col_major>::set_local)
    .def("set_global", &wrap_distributed_matrix<matrix_col_major>::set_global)
    .def("update_local", &wrap_distributed_matrix<matrix_col_major>::update_local)
    .def("update_global", &wrap_distributed_matrix<matrix_col_major>::update_global)
    .def("set_zeros", &wrap_distributed_matrix<matrix_col_major>::set_zeros)
    .def_property_readonly("length_array", &wrap_distributed_matrix<matrix_col_major>::get_length_array)
    .def_property_readonly("lld", &wrap_distributed_matrix<matrix_col_major>::get_lld)
    //.def("generate", &wrap_distributed_matrix<matrix_col_major>::generate)
    .def("print", &wrap_distributed_matrix<matrix_col_major>::print)
    .def("get_map", &wrap_distributed_matrix<matrix_col_major>::get_map)
    .def_property_readonly("major", &wrap_distributed_matrix<matrix_col_major>::get_major_string)
    .def("set_ndarray", &wrap_distributed_matrix<matrix_col_major>::set_ndarray)
    .def("get_ndarray", &wrap_distributed_matrix<matrix_col_major>::get_ndarray, py::return_value_policy::reference_internal)
    .def_property("ndarray", &wrap_distributed_matrix<matrix_col_major>::get_ndarray, &wrap_distributed_matrix<matrix_col_major>::set_ndarray)
    .def_property_readonly("map", &wrap_distributed_matrix<matrix_col_major>::get_map)
    .def_property_readonly("major", &wrap_distributed_matrix<matrix_col_major>::get_major_string);

  m.def("product", py::overload_cast<double, wrap_distributed_matrix<matrix_col_major> const&, bool, wrap_distributed_matrix<matrix_col_major> const&, bool, double, wrap_distributed_matrix<matrix_col_major>&>(&pyrokko_product), py::arg("alpha"), py::arg("matA"), py::arg("transA"), py::arg("matB"), py::arg("transB"), py::arg("beta"), py::arg("matC"));

  py::class_<wrap_parallel_dense_ev>(m, "parallel_dense_ev")
    .def(py::init<std::string const&>())
    .def(py::init<>())
    .def("initialize", &wrap_parallel_dense_ev::initialize)
    .def("finalize", &wrap_parallel_dense_ev::finalize)
    .def("default_mapping", &wrap_parallel_dense_ev::default_mapping)
    .def("diagonalize", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&, Eigen::RefVec<double>&, wrap_parameters const&>(&wrap_parallel_dense_ev::diagonalize<Eigen::RefVec<double>>),
         py::arg("mat"), py::arg("eigvals"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&, Eigen::RefVec<double>&, wrap_distributed_matrix<matrix_col_major>&, wrap_parameters const&>(&wrap_parallel_dense_ev::diagonalize<Eigen::RefVec<double>>),
         py::arg("mat"), py::arg("eigvals"), py::arg("eigvecs"), py::arg("params") = wrap_parameters())
    .def_property_readonly_static("solvers", &wrap_parallel_dense_ev::solvers)
    .def_property_readonly_static("default_solver", &wrap_parallel_dense_ev::default_solver);


  py::class_<wrap_minij_matrix>(m, "minij_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_minij_matrix::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_minij_matrix::generate_col))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&>(&wrap_minij_matrix::generate<matrix_col_major>))
    .def_static("eigenvalue", &minij_matrix::eigenvalue);
  
  py::class_<wrap_frank_matrix>(m, "frank_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_frank_matrix::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_frank_matrix::generate_col))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&>(&wrap_frank_matrix::generate<matrix_col_major>))
    .def_static("eigenvalue", &frank_matrix::eigenvalue);

  py::class_<wrap_laplacian_matrix>(m, "laplacian_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_laplacian_matrix::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_laplacian_matrix::generate_col))
    .def_static("eigenvalue", &laplacian_matrix::eigenvalue);

  py::class_<wrap_helmert_matrix>(m, "helmert_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_helmert_matrix::generate_eigen<Eigen::RowMajor>))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_helmert_matrix::generate_eigen<Eigen::ColMajor>))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&>(&wrap_helmert_matrix::generate<matrix_col_major>))
    .def_static("generate_for_given_eigenvalues", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>>(&wrap_helmert_matrix::generate_for_given_eigenvalues_eigen<Eigen::RowMajor>))
    .def_static("generate_for_given_eigenvalues", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>>(&wrap_helmert_matrix::generate_for_given_eigenvalues_eigen<Eigen::ColMajor>))
    .def_static("generate_for_given_eigenvalues", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&,Eigen::Ref<Eigen::VectorXd>>(&wrap_helmert_matrix::generate_for_given_eigenvalues<matrix_col_major>));

  py::class_<wrap_matrix012>(m, "matrix012")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_matrix012::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_matrix012::generate_col))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix<matrix_col_major>&>(&wrap_matrix012::generate<matrix_col_major>));

  // collective MPI communication
  m.def("gather", py::overload_cast<wrap_distributed_matrix<matrix_row_major> const&, Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>, int>(&pyrokko_gather));
  m.def("gather", py::overload_cast<wrap_distributed_matrix<matrix_col_major> const&, Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>, int>(&pyrokko_gather));

  m.def("scatter", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>, wrap_distributed_matrix<matrix_row_major>&, int>(&pyrokko_scatter));
  m.def("scatter", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>, wrap_distributed_matrix<matrix_col_major>&, int>(&pyrokko_scatter));

  py::class_<distributed_mfree>(m, "distributed_mfree_base");

  py::class_<distributed_mfree_default, distributed_mfree>(m, "distributed_mfree_default");

  py::class_<distributed_mfree_inherit, wrap_distributed_mfree_inherit, distributed_mfree_default, distributed_mfree>(m, "distributed_mfree_inherit")
    .def(py::init<int>())
    .def(py::init<int, pybind11::handle const&>())
    .def("multiply", py::overload_cast<const Eigen::VectorXd&,Eigen::Ref<Eigen::VectorXd>>(&distributed_mfree_inherit::multiply, py::const_))
    .def_property_readonly("dim", &distributed_mfree_default::get_dim)
    .def_property_readonly("num_local_rows", &distributed_mfree_default::get_num_local_rows)
    .def_property_readonly("start_row", &distributed_mfree_default::start_row)
    .def_property_readonly("end_row", &distributed_mfree_default::end_row);

  py::class_<wrap_distributed_mfree>(m, "distributed_mfree")
    .def(py::init<std::function<void(ConstMapVec,MapVec)>, int>())
    .def(py::init<std::function<void(ConstMapVec,MapVec)>, int, pybind11::handle const&>())
    .def(py::init<std::function<void(ConstMapVec,MapVec)>, wrap_mapping_1d const&>())
    .def_property_readonly("dim", &distributed_mfree_default::get_dim)
    .def_property_readonly("num_local_rows", &distributed_mfree_default::get_num_local_rows)
    .def_property_readonly("start_row", &distributed_mfree_default::start_row)
    .def_property_readonly("end_row", &distributed_mfree_default::end_row);

  // sparse
  py::class_<wrap_mapping_1d>(m, "mapping_1d")
    .def(py::init<int, pybind11::handle const&, std::string const&>())
    .def_property_readonly("dim", &mapping_1d::get_dim)
    .def_property_readonly("num_local_rows", &mapping_1d::num_local_rows)
    .def_property_readonly("start_row", &mapping_1d::start_row)
    .def_property_readonly("end_row", &mapping_1d::end_row)
    .def_property_readonly("get_solver_name", &mapping_1d::get_solver_name);

  py::class_<wrap_distributed_crs_matrix>(m, "distributed_crs_matrix")
    .def(py::init<wrap_mapping_1d, int>())
    .def("insert", py::overload_cast<int, std::vector<int> const&, std::vector<double> const&>(&wrap_distributed_crs_matrix::insert, py::const_))
    .def("complete", &distributed_crs_matrix::complete)
    .def_property_readonly("dim", &distributed_crs_matrix::get_dim)
    .def_property_readonly("num_local_rows", &distributed_crs_matrix::get_num_local_rows)
    .def_property_readonly("start_row", &distributed_crs_matrix::start_row)
    .def_property_readonly("end_row", &distributed_crs_matrix::end_row)
    .def_property_readonly("nnz", &distributed_crs_matrix::get_nnz)
    .def("print", &distributed_crs_matrix::print)
    .def("output_matrix_market", &distributed_crs_matrix::output_matrix_market)
    .def_property_readonly("solver_name", &distributed_crs_matrix::get_solver_name);

  py::class_<wrap_parallel_sparse_ev>(m, "parallel_sparse_ev")
    .def(py::init<std::string const&>())
    .def(py::init<>())
    .def("initialize", &wrap_parallel_sparse_ev::initialize)
    .def("finalize", &parallel_sparse_ev::finalize)
    .def("default_mapping", &wrap_parallel_sparse_ev::default_mapping)
    .def("eigenvalue", &parallel_sparse_ev::eigenvalue)
    .def("eigenvector", &wrap_parallel_sparse_ev::python_eigenvector)
    .def_property_readonly("num_conv", &parallel_sparse_ev::get_num_conv)
    .def("diagonalize", py::overload_cast<wrap_distributed_crs_matrix const&, wrap_parameters const&>(&wrap_parallel_sparse_ev::diagonalize),
         py::arg("mat"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<distributed_mfree const&, wrap_parameters const&>(&wrap_parallel_sparse_ev::diagonalize),
         py::arg("mat"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<wrap_distributed_mfree const&, wrap_parameters const&>(&wrap_parallel_sparse_ev::diagonalize),
         py::arg("mat"), py::arg("params") = wrap_parameters())
    .def_property_readonly_static("solvers", &parallel_sparse_ev::solvers)
    .def_property_readonly_static("default_solver", &parallel_sparse_ev::default_solver);
  
  // Eigen3 matrix
  m.def("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,std::function<double(int, int)>const&>(&pyrokko_generate<Eigen::ColMajor>));
  m.def("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,std::function<double(int, int)>const&>(&pyrokko_generate<Eigen::RowMajor>));

  // utility functions
  m.def("sort_eigenpairs", py::overload_cast<Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,bool>(&pyrokko_sort_eigenpairs<Eigen::RowMajor>),
        py::arg("eigval"), py::arg("eigvec"), py::arg("eigval_sorted"), py::arg("eigvec_sorted"), py::arg("ascending") = true);
  
  m.def("sort_eigenpairs", py::overload_cast<Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,bool>(&pyrokko_sort_eigenpairs<Eigen::ColMajor>),
        py::arg("eigval"), py::arg("eigvec"), py::arg("eigval_sorted"), py::arg("eigvec_sorted"), py::arg("ascending") = true);

  m.def("split_solver_name", &wrap_split_solver_name);
}

} // end namespace rokko
