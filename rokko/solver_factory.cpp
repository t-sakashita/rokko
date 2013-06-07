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
