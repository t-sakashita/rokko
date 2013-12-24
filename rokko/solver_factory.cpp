/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
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

solver_factory::solver_pointer_type solver_factory::make_solver(std::string const& name) {
  return instance()->make_creator(name)->create();
}

solver_factory::creator_pointer_type solver_factory::make_creator(std::string const& name) const {
  creator_map_type::const_iterator itr = creators_.find(name);
  if (itr == creators_.end() || itr->second == 0) {
    std::cerr << "Error: unknown solver: \"" << name << "\" (registered solvers: ";
    for (creator_map_type::const_iterator itr = creators_.begin(); itr != creators_.end(); ++itr) {
      if (itr != creators_.begin()) std::cerr << ", ";
      std::cerr << "\"" << itr->first << "\"";
    }
    std::cerr << ")\n";
    boost::throw_exception(std::runtime_error("solver_factory::make_creator()"));
  }
  return itr->second;
}

std::vector<std::string> solver_factory::solver_names() const {
    std::vector<std::string> retvec;
    for (creator_map_type::const_iterator it = creators_.begin();
         it != creators_.end(); ++it) {
        retvec.push_back(it->first);
    }
    return retvec;
}

solver_factory* solver_factory::instance() {
  if (!instance_) instance_ = new solver_factory;
  return instance_;
}

solver_factory* solver_factory::instance_ = 0;

} // end namespace rokko

#ifdef ROKKO_BUILE_SHARED_LIBS
# include <rokko/scalapack/core.hpp>
# include <rokko/eigen_exa/core.hpp>
# include <rokko/elpa/core.hpp>
# include <rokko/elemental/core.hpp>
  ROKKO_REGISTER_SOLVER(rokko::scalapack::solver<rokko::scalapack::pdsyev>, "scalapack")
  ROKKO_REGISTER_SOLVER(rokko::eigen_exa::solver, "eigen_exa")
  ROKKO_REGISTER_SOLVER(rokko::elpa::solver, "elpa")
  ROKKO_REGISTER_SOLVER(rokko::elemental::solver, "elemental")
#endif
