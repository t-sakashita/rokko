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

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>
#include <boost/throw_exception.hpp>

namespace rokko {

class parameters {
private:
  typedef std::map<std::string, boost::any> map_type;
  typedef map_type::mapped_type mapped_type;
public:
  typedef map_type::key_type key_type;
  typedef map_type::value_type value_type;

  parameters() {};
  ~parameters() {};
  
  void clear() { map_.clear(); }
  void clear(key_type const& key) { map_.erase(map_.find(key)); }
  // set "value" for "key"
  template<typename T>
  void set(key_type const& key, T const& value) {
    map_[key] = value;
  }
  void set(key_type const& key, const char *value) {
    map_[key] = std::string(value);
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
      BOOST_THROW_EXCEPTION(std::invalid_argument("parameters::get() : type given as template parameter is not correct."));
    }
    return boost::any_cast<T>(map_.find(key)->second);
  }
  const std::type_info& type(key_type const& key) const { return map_.find(key)->second.type(); }
  // return parameter as string
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
      BOOST_THROW_EXCEPTION(std::invalid_argument("parameters::get_string() : value type given as template parameter must be char*, string, int, or double.")); 
    }
  }
  bool get_bool(key_type const& key) const {
    if (defined(key)) {
      if (type(key) == typeid(bool)) {
	return get<bool>(key);
      }
      else {
	BOOST_THROW_EXCEPTION(std::invalid_argument("parameters::get_bool() : the key \"" + key + "\" is not bool type")); 
      }
    }
    return false;
  }
  std::map<std::string, boost::any>& get_map() { return map_; }
  
private:
  std::map<std::string, boost::any> map_;
};


//void get_2string(rokko::parameters const& params, std::string const& key, char keyword_char, std::string const& keyword_str) {
//  if (params.type(key) == typeid(const char*)) {
//    letter = *(params.get<const char*>(key));
//  } else if (params.type(key) == typeid(const char*)) {
//    letter = params.get<char>(key);    
//  }
//}

template <typename T>
bool get_key(rokko::parameters const& params, std::string const& key, T& value) {
  if (params.defined(key)) {
    if (params.type(key) == typeid(T)) {
      value = params.get<T>(key);
      return true;
    }
  }
  return false;
}


} // namespace rokko

#endif // ROKKO_PARAMETERS_HPP
