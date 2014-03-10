/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*               2013-2013    Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include "solver_factory.hpp"
#include <stdexcept>

namespace rokko {

solver_factory::serial_dense_solver_pointer_type solver_factory::make_serial_dense_solver(
  std::string const& name) {
  solver_factory* factory = solver_factory::instance();
  if (name == "") {
    if (instance()->default_sd_solver_ != "") {
      return factory->make_sd_creator(factory->default_sd_solver_)->create();
    } else {
      std::cerr << "Error: default solver is not defined\n";
      boost::throw_exception(std::runtime_error("solver_factory::make_serial_dense_solver()"));
    }
  } else {
    return factory->make_sd_creator(name)->create();
  }
}

solver_factory::serial_dense_solver_pointer_type solver_factory::make_serial_dense_solver() {
  return solver_factory::make_serial_dense_solver("");
}

solver_factory::sd_creator_pointer_type
solver_factory::make_sd_creator(std::string const& name) const {
  sd_creator_map_type::const_iterator itr = sd_creators_.find(name);
  if (itr == sd_creators_.end() || itr->second == 0) {
    std::cerr << "Error: unknown solver: \"" << name << "\" (registered solvers: ";
    for (sd_creator_map_type::const_iterator itr = sd_creators_.begin();
         itr != sd_creators_.end(); ++itr) {
      if (itr != sd_creators_.begin()) std::cerr << ", ";
      std::cerr << "\"" << itr->first << "\"";
    }
    std::cerr << ")\n";
    boost::throw_exception(std::runtime_error("solver_factory::make_sd_creator()"));
  }
  return itr->second;
}

solver_factory::parallel_dense_solver_pointer_type solver_factory::make_parallel_dense_solver(
  std::string const& name) {
  solver_factory* factory = solver_factory::instance();
  if (name == "") {
    if (instance()->default_pd_solver_ != "") {
      return factory->make_pd_creator(factory->default_pd_solver_)->create();
    } else {
      std::cerr << "Error: default solver is not defined\n";
      boost::throw_exception(std::runtime_error("solver_factory::make_parallel_dense_solver()"));
    }
  } else {
    return factory->make_pd_creator(name)->create();
  }
}

solver_factory::parallel_dense_solver_pointer_type solver_factory::make_parallel_dense_solver() {
  return solver_factory::make_parallel_dense_solver("");
}

solver_factory::pd_creator_pointer_type
solver_factory::make_pd_creator(std::string const& name) const {
  pd_creator_map_type::const_iterator itr = pd_creators_.find(name);
  if (itr == pd_creators_.end() || itr->second == 0) {
    std::cerr << "Error: unknown solver: \"" << name << "\" (registered solvers: ";
    for (pd_creator_map_type::const_iterator itr = pd_creators_.begin();
         itr != pd_creators_.end(); ++itr) {
      if (itr != pd_creators_.begin()) std::cerr << ", ";
      std::cerr << "\"" << itr->first << "\"";
    }
    std::cerr << ")\n";
    boost::throw_exception(std::runtime_error("solver_factory::make_pd_creator()"));
  }
  return itr->second;
}

std::vector<std::string> solver_factory::serial_dense_solver_names() {
  solver_factory* factory = solver_factory::instance();
  std::vector<std::string> retvec;
  for (sd_creator_map_type::const_iterator it = factory->sd_creators_.begin();
       it != factory->sd_creators_.end(); ++it) {
    retvec.push_back(it->first);
  }
  return retvec;
}

std::string solver_factory::default_serial_dense_solver_name() {
  return instance()->default_sd_solver_;
}

std::vector<std::string> solver_factory::parallel_dense_solver_names() {
  solver_factory* factory = solver_factory::instance();
  std::vector<std::string> retvec;
  for (pd_creator_map_type::const_iterator it = factory->pd_creators_.begin();
       it != factory->pd_creators_.end(); ++it) {
    retvec.push_back(it->first);
  }
  return retvec;
}

std::string solver_factory::default_parallel_dense_solver_name() {
  return instance()->default_pd_solver_;
}

solver_factory* solver_factory::instance() {
  if (!instance_) instance_ = new solver_factory;
  return instance_;
}

solver_factory* solver_factory::instance_ = 0;

} // end namespace rokko

#ifndef ROKKO_BUILD_SHARED_LIBS
#ifdef BUILD_SCALAPACK
# include <rokko/scalapack/core.hpp>
  ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(rokko::scalapack::solver<rokko::scalapack::pdsyev>, "scalapack", 20)
#endif
#ifdef BUILD_EIGENEXA
# include <rokko/eigen_exa/core.hpp>
  ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(rokko::eigen_exa::solver, "eigen_exa", 40)
#endif
#ifdef BUILD_ELPA
# include <rokko/elpa/core.hpp>
  ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(rokko::elpa::solver, "elpa", 30)
#endif
#ifdef BUILD_ELEMENTAL
# include <rokko/elemental/core.hpp>
  ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(rokko::elemental::solver, "elemental", 10)
#endif
#endif
