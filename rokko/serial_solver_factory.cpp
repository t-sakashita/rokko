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

#include "serial_solver_factory.hpp"
#include <stdexcept>

namespace rokko {

serial_solver_factory::solver_pointer_type serial_solver_factory::make_solver(std::string const& name) {
  return instance()->make_creator(name)->create();
}

serial_solver_factory::creator_pointer_type serial_solver_factory::make_creator(std::string const& name) const {
  creator_map_type::const_iterator itr = creators_.find(name);
  if (itr == creators_.end() || itr->second == 0) {
    std::cerr << "Error: unknown solver: \"" << name << "\" (registered solvers: ";
    for (creator_map_type::const_iterator itr = creators_.begin(); itr != creators_.end(); ++itr) {
      if (itr != creators_.begin()) std::cerr << ", ";
      std::cerr << "\"" << itr->first << "\"";
    }
    std::cerr << ")\n";
    boost::throw_exception(std::runtime_error("serial_solver_factory::make_creator()"));
  }
  return itr->second;
}

std::vector<std::string> serial_solver_factory::solver_names() {
  serial_solver_factory* factory = serial_solver_factory::instance();
  std::vector<std::string> retvec;
  for (creator_map_type::const_iterator it = factory->creators_.begin();
       it != factory->creators_.end(); ++it) {
    retvec.push_back(it->first);
  }
  return retvec;
}

serial_solver_factory* serial_solver_factory::instance() {
  if (!instance_) instance_ = new serial_solver_factory;
  return instance_;
}

serial_solver_factory* serial_solver_factory::instance_ = 0;

} // end namespace rokko
