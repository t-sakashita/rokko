/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PARAMETERS_HPP
#define ROKKO_PARAMETERS_HPP

#include <list>
#include <map>
#include <string>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>

namespace rokko {

class parameters {
private:
  typedef std::map<std::string, boost::any> map_type;
  typedef map_type::value_type value_type;
  typedef map_type::mapped_type mapped_type;
public:
  typedef map_type::key_type key_type;
  void clear() { map_.clear(); }
  void clear(key_type const& key) { map_.erase(map_.find(key)); }
  // set "value" for "key"
  template<typename T>
  void set(key_type const& key, T const& value) {
    map_[key] = value;
  }
  // returns if parameter with "key" is defined or not
  bool defined(key_type const& key) const {return (map_.find(key) != map_.end()); }
  // returns list of keys
  std::list<key_type> keys() const {
    std::list<key_type> keys;
    BOOST_FOREACH(value_type const& p, map_) { keys.push_back(p.first); }
    return keys;
  }
  // returns value of parameter in type T
  template<typename T>
  T get(key_type const& key) const {
    if (type(key) != typeid(T)) {
      throw "error: type";
    }
    return boost::any_cast<T>(map_.find(key)->second);
  }
  const std::type_info& type(key_type const& key) const { return map_.find(key)->second.type(); }
  std::string get_string(key_type const& key) const {
    if (type(key) == typeid(std::string)) {
      return get<std::string>(key);
    }
    if (type(key) == typeid(const char*)) {
      return std::string(get<const char*>(key));
    }
    if (type(key) == typeid(int)) {
      return boost::lexical_cast<std::string>(get<int>(key));
    }
    if (type(key) == typeid(double)) {
      return boost::lexical_cast<std::string>(get<double>(key));
    }
    if (type(key) == typeid(bool)) {
      return boost::lexical_cast<std::string>(get<bool>(key));
    }
    if (type(key) == typeid(char)) {
      return boost::lexical_cast<std::string>(get<char>(key));
    }
    else {
      throw "error: set_string only accepts charatcters, string, int or double as an argument";      
    }
  }
private:
  std::map<std::string, boost::any> map_;
};

} // namespace rokko

#endif // ROKKO_PARAMETERS_HPP
