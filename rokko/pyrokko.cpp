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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <rokko/eigen3.hpp>
#include <rokko/pyrokko_serial_dense_ev.hpp>

#include <rokko/pyrokko_parameters.hpp>
#include <rokko/pyrokko_grid.hpp>
#include <rokko/pyrokko_mapping_bc.hpp>
#include <rokko/pyrokko_distributed_matrix.hpp>
#include <rokko/pyrokko_parallel_dense_ev.hpp>
#include <rokko/pyrokko_collective.hpp>

#include <rokko/pyrokko_localized_matrix.hpp>

#include <rokko/utility/pyrokko_minij_matrix.hpp>
#include <rokko/utility/pyrokko_frank_matrix.hpp>
#include <rokko/utility/pyrokko_laplacian_matrix.hpp>
#include <rokko/utility/pyrokko_matrix012.hpp>

#include <rokko/pyrokko_parallel_sparse_ev.hpp>
#include <rokko/distributed_crs_matrix.hpp>
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
  
  py::class_<wrap_localized_matrix>(m, "localized_matrix")
    .def(py::init<matrix_major_enum>(), py::arg("major") = col)
    .def(py::init<int, int, matrix_major_enum>(), py::arg("rows"), py::arg("cols"), py::arg("major") = col)
    .def("generate", &wrap_localized_matrix::generate)
    .def_property_readonly("major", &wrap_localized_matrix::get_major_string);


  py::class_<wrap_serial_dense_ev>(m, "serial_dense_ev")
    .def(py::init<std::string const&>())
    .def(py::init<>())
    //.def("initialize", py::overload_cast<int&, char**&>(&serial_dense_ev::initialize))
    .def("finalize", &serial_dense_ev::finalize)
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<Eigen::ColMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("eigvec"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::VectorXd>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<Eigen::RowMajor>),
         py::arg("mat"), py::arg("eigval"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::VectorXd>,wrap_parameters const&>(&wrap_serial_dense_ev::diagonalize<Eigen::ColMajor>),
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
    .def(py::init<pybind11::handle const&, grid_row_major_t, int>())
    .def(py::init<pybind11::handle const&, grid_col_major_t, int>())
    .def("get_comm", &wrap_grid::get_comm)
    .def_property_readonly("nprocs", &wrap_grid::get_nprocs)
    .def_property_readonly("nprow", &wrap_grid::get_nprow)
    .def_property_readonly("npcol", &wrap_grid::get_npcol)
    .def_property_readonly("myrank", &wrap_grid::get_myrank)
    .def_property_readonly("myrow", &wrap_grid::get_myrow)
    .def_property_readonly("mycol", &wrap_grid::get_mycol)
    .def_property_readonly("shape", &wrap_grid::get_shape)
    .def_property_readonly("mine", &wrap_grid::get_mine)
    .def_property_readonly("major", &wrap_grid::get_major_string);

  py::class_<wrap_mapping_bc>(m, "mapping_bc")
    .def(py::init<matrix_major_enum>(), py::arg("major") = col)
    .def(py::init<int, int, wrap_grid, matrix_major_enum>(), py::arg("global_dim"), py::arg("block_size"), py::arg("grid"), py::arg("major") = col)
    .def(py::init<int, int, int, wrap_grid, matrix_major_enum>(), py::arg("global_dim"), py::arg("block_size"), py::arg("lld"), py::arg("grid"), py::arg("major") = col)
    .def("get_mb", &wrap_mapping_bc::get_mb)
    .def("get_nb", &wrap_mapping_bc::get_nb)
    .def_property_readonly("block_shape", &wrap_mapping_bc::get_block_shape)
    .def_property_readonly("m_global", &wrap_mapping_bc::get_m_global)
    .def_property_readonly("n_global", &wrap_mapping_bc::get_n_global)
    .def_property_readonly("global_shape", &wrap_mapping_bc::get_global_shape)
    .def_property_readonly("m_local", &wrap_mapping_bc::get_m_local)
    .def_property_readonly("n_local", &wrap_mapping_bc::get_n_local)
    .def_property_readonly("local_shape", &wrap_mapping_bc::get_local_shape)
    .def("translate_l2g_row", &wrap_mapping_bc::translate_l2g_row)
    .def("translate_l2g_col", &wrap_mapping_bc::translate_l2g_col)
    .def("translate_l2g", &wrap_mapping_bc::translate_l2g) // tuple
    .def("translate_g2l_row", &wrap_mapping_bc::translate_g2l_row)
    .def("translate_g2l_col", &wrap_mapping_bc::translate_g2l_col)
    .def("translate_g2l", &wrap_mapping_bc::translate_g2l) // tuple
    .def_property_readonly("is_gindex_myrow", &wrap_mapping_bc::is_gindex_myrow)
    .def_property_readonly("is_gindex_mycol", &wrap_mapping_bc::is_gindex_mycol)
    .def_property_readonly("is_gindex", &wrap_mapping_bc::is_gindex)
    .def_property_readonly("major", &wrap_mapping_bc::get_major_string);
  
  py::class_<wrap_distributed_matrix>(m, "distributed_matrix")
    .def(py::init<wrap_mapping_bc>())
    .def_property_readonly("mb", &wrap_distributed_matrix::get_mb)
    .def_property_readonly("nb", &wrap_distributed_matrix::get_nb)
    .def_property_readonly("block_shape", &wrap_distributed_matrix::get_block_shape)
    .def_property_readonly("m_global", &wrap_distributed_matrix::get_m_global)
    .def_property_readonly("n_global", &wrap_distributed_matrix::get_n_global)
    .def_property_readonly("m_local", &wrap_distributed_matrix::get_m_local)
    .def_property_readonly("n_local", &wrap_distributed_matrix::get_n_local)
    .def_property_readonly("global_shape", &wrap_distributed_matrix::get_global_shape)
    .def_property_readonly("local_shape", &wrap_distributed_matrix::get_local_shape)
    .def("translate_l2g_row", &wrap_distributed_matrix::translate_l2g_row)
    .def("translate_l2g_col", &wrap_distributed_matrix::translate_l2g_col)
    .def("translate_l2g", &wrap_distributed_matrix::translate_l2g) // tuple
    .def("translate_g2l_row", &wrap_distributed_matrix::translate_g2l_row)
    .def("translate_g2l_col", &wrap_distributed_matrix::translate_g2l_col)
    .def("translate_g2l", &wrap_distributed_matrix::translate_g2l) // tuple
    .def_property_readonly("is_gindex_myrow", &wrap_distributed_matrix::is_gindex_myrow)
    .def_property_readonly("is_gindex_mycol", &wrap_distributed_matrix::is_gindex_mycol)
    .def_property_readonly("is_gindex", &wrap_distributed_matrix::is_gindex)
    .def("set_local", &wrap_distributed_matrix::set_local)
    .def("set_global", &wrap_distributed_matrix::set_global)
    .def("set_zeros", &wrap_distributed_matrix::set_zeros)
    .def_property_readonly("length_array", &wrap_distributed_matrix::get_length_array)
    .def_property_readonly("lld", &wrap_distributed_matrix::get_lld)
    .def("generate", &wrap_distributed_matrix::generate)
    .def("print", &wrap_distributed_matrix::print)
    .def("get_map", &wrap_distributed_matrix::get_map)
    .def_property_readonly("major", &wrap_distributed_matrix::get_major_string)
    .def("set_ndarray", &wrap_distributed_matrix::set_ndarray)
    .def("get_ndarray", &wrap_distributed_matrix::get_ndarray, py::return_value_policy::reference_internal)
    .def_property("ndarray", &wrap_distributed_matrix::get_ndarray, &wrap_distributed_matrix::set_ndarray)
    .def_property_readonly("major", &wrap_distributed_matrix::get_major_string);


  py::class_<wrap_parallel_dense_ev>(m, "parallel_dense_ev")
    .def(py::init<std::string const&>())
    .def(py::init<>())
    .def("initialize", &wrap_parallel_dense_ev::initialize)
    .def("finalize", &wrap_parallel_dense_ev::finalize)
    .def("default_mapping", &wrap_parallel_dense_ev::default_mapping)
    .def("diagonalize", py::overload_cast<wrap_distributed_matrix&, Eigen::RefVec<double>&, wrap_parameters const&>(&wrap_parallel_dense_ev::diagonalize<Eigen::RefVec<double>>),
         py::arg("mat"), py::arg("eigvals"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<wrap_distributed_matrix&, Eigen::RefVec<double>&, wrap_distributed_matrix&, wrap_parameters const&>(&wrap_parallel_dense_ev::diagonalize<Eigen::RefVec<double>>),
         py::arg("mat"), py::arg("eigvals"), py::arg("eigvecs"), py::arg("params") = wrap_parameters())
    .def_property_readonly_static("solvers", &wrap_parallel_dense_ev::solvers)
    .def_property_readonly_static("default_solver", &wrap_parallel_dense_ev::default_solver);


  py::class_<wrap_minij_matrix>(m, "minij_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_minij_matrix::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_minij_matrix::generate_col))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix&>(&wrap_minij_matrix::generate))
    .def_static("eigenvalue", &minij_matrix::eigenvalue);
  
  py::class_<wrap_frank_matrix>(m, "frank_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_frank_matrix::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_frank_matrix::generate_col))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix&>(&wrap_frank_matrix::generate))
    .def_static("eigenvalue", &frank_matrix::eigenvalue);

  py::class_<wrap_laplacian_matrix>(m, "laplacian_matrix")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_laplacian_matrix::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_laplacian_matrix::generate_col))
    .def_static("eigenvalue", &laplacian_matrix::eigenvalue);

  py::class_<wrap_matrix012>(m, "matrix012")
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>>(&wrap_matrix012::generate_row))
    .def_static("generate", py::overload_cast<Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>>(&wrap_matrix012::generate_col))
    .def_static("generate", py::overload_cast<wrap_distributed_matrix&>(&wrap_matrix012::generate));

  // collective MPI communication
  m.def("gather", &pyrokko_gather);
  m.def("scatter", &pyrokko_scatter);

  py::class_<distributed_mfree>(m, "distributed_mfree_base");

  py::class_<wrap_distributed_mfree>(m, "distributed_mfree")
    .def(py::init<std::function<void(ConstMapVec,MapVec)>, int, int>());
  
  // sparse
  py::class_<wrap_parallel_sparse_ev>(m, "parallel_sparse_ev")
    .def(py::init<std::string const&>())
    .def(py::init<>())
    .def("initialize", &wrap_parallel_sparse_ev::initialize)
    .def("finalize", &parallel_sparse_ev::finalize)
    .def("eigenvalue", &parallel_sparse_ev::eigenvalue)
    .def("eigenvector", &wrap_parallel_sparse_ev::python_eigenvector)
    .def_property_readonly("num_conv", &parallel_sparse_ev::num_conv)
    .def("diagonalize", py::overload_cast<distributed_crs_matrix&, wrap_parameters const&>(&wrap_parallel_sparse_ev::diagonalize),
         py::arg("mat"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<distributed_mfree*, wrap_parameters const&>(&wrap_parallel_sparse_ev::diagonalize),
         py::arg("mat"), py::arg("params") = wrap_parameters())
    .def("diagonalize", py::overload_cast<wrap_distributed_mfree&, wrap_parameters const&>(&wrap_parallel_sparse_ev::diagonalize),
         py::arg("mat"), py::arg("params") = wrap_parameters())
    .def_property_readonly_static("solvers", &parallel_sparse_ev::solvers)
    .def_property_readonly_static("default_solver", &parallel_sparse_ev::default_solver);
  
  py::class_<distributed_crs_matrix>(m, "distributed_crs_matrix")
    .def(py::init<int, int, wrap_parallel_sparse_ev&>())
    .def(py::init<int, int, int, wrap_parallel_sparse_ev&>())
    .def("insert", py::overload_cast<int, std::vector<int> const&, std::vector<double> const&>(&distributed_crs_matrix::insert))
    .def("complete", &distributed_crs_matrix::complete)
    .def_property_readonly("dim", &distributed_crs_matrix::get_dim)
    .def_property_readonly("num_local_rows", &distributed_crs_matrix::num_local_rows)
    .def_property_readonly("start_row", &distributed_crs_matrix::start_row)
    .def_property_readonly("end_row", &distributed_crs_matrix::end_row)
    .def_property_readonly("nnz", &distributed_crs_matrix::get_nnz)
    .def("print", &distributed_crs_matrix::print)
    .def("output_matrix_market", &distributed_crs_matrix::output_matrix_market)
    .def_property_readonly("solver_name", &distributed_crs_matrix::get_solver_name);

  // utility functions
  m.def("sort_eigenpairs", py::overload_cast<Eigen::Ref<Eigen::Vector<double, Eigen::Dynamic>>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,Eigen::Ref<Eigen::Vector<double, Eigen::Dynamic>>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>,bool>(&pyrokko_sort_eigenpairs<Eigen::RowMajor>),
        py::arg("eigval"), py::arg("eigvec"), py::arg("eigval_sorted"), py::arg("eigvec_sorted"), py::arg("ascending") = true);
  
  m.def("sort_eigenpairs", py::overload_cast<Eigen::Ref<Eigen::Vector<double, Eigen::Dynamic>>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,Eigen::Ref<Eigen::Vector<double, Eigen::Dynamic>>,Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>,bool>(&pyrokko_sort_eigenpairs<Eigen::ColMajor>),
        py::arg("eigval"), py::arg("eigvec"), py::arg("eigval_sorted"), py::arg("eigvec_sorted"), py::arg("ascending") = true);

  m.def("split_solver_name", &wrap_split_solver_name);
}

} // end namespace rokko
