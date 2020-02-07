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

#ifndef PYROKKO_PARAMETERS_HPP
#define PYROKKO_PARAMETERS_HPP

#include <rokko/parameters.hpp>

namespace rokko {

namespace py = pybind11;

class wrap_parameters : public rokko::parameters {
public:
  wrap_parameters() = default;
  wrap_parameters(rokko::parameters const& params_in) : parameters(params_in) {}
  py::object python_get(std::string const& key) const {
    if (type(key) == typeid(int)) {
      return py::cast(get<int>(key));
    } else if (type(key) == typeid(double)) {
      return py::cast(get<double>(key));
    } else if (type(key) == typeid(std::string)) {
      return py::cast(get<std::string>(key));
    } else if (type(key) == typeid(bool)) {
      return py::cast(get<bool>(key));
    } else if (type(key) == typeid(char)) {
      return py::cast(get<char>(key));
    }
    throw std::invalid_argument("wrap_parameters::python_get() : value type given as template parameter must be char*, string, int, or double.");
  }
  py::list python_keys() const {
    py::list py_list;
    for(auto const& p : get_map()) {
      py_list.append(p.first);
    }
    return py_list;
  }
  py::dict dict() const {
    py::dict dict;
    for(auto const& p : get_map()) {
      dict[p.first.c_str()] = python_get(p.first);
    }
    return dict;
  }
};

} // end namespace rokko

#endif // PYROKKO_PARAMETERS_HPP
