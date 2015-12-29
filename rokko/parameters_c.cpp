/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parameters.h>
#include <rokko/parameters.hpp>

void rokko_parameters_construct(struct rokko_parameters* params) {
  params->ptr = new rokko::parameters();
}

void rokko_parameters_destruct(struct rokko_parameters* params) {
  rokko::parameters* ptr = static_cast<rokko::parameters*>(params->ptr);
  delete ptr;
  params->ptr = 0;
}

void rokko_parameters_clear_all(struct rokko_parameters params) {
  static_cast<rokko::parameters*>(params.ptr)->clear();
}

void rokko_parameters_clear(struct rokko_parameters params, const char* key) {
  static_cast<rokko::parameters*>(params.ptr)->clear(key);
}

void rokko_parameters_set_int(rokko_parameters params, const char* key, int value) {
  static_cast<rokko::parameters*>(params.ptr)->set<int>(key, value);
}

void rokko_parameters_set_double(rokko_parameters params, const char* key, double value) { 
  static_cast<rokko::parameters*>(params.ptr)->set<double>(key, value);
}

void rokko_parameters_set_true(rokko_parameters params, const char* key) {
  static_cast<rokko::parameters*>(params.ptr)->set<bool>(key, true);
}

void rokko_parameters_set_false(rokko_parameters params, const char* key) {
  static_cast<rokko::parameters*>(params.ptr)->set<bool>(key, false);
}

void rokko_parameters_set_char(rokko_parameters params, const char* key, char value) { 
  static_cast<rokko::parameters*>(params.ptr)->set<char>(key, value);
}

void rokko_parameters_set_string(struct rokko_parameters params, const char* key, const char* value) {
  static_cast<rokko::parameters*>(params.ptr)->set<std::string>(key, std::string(value));
}

int rokko_parameters_get_int(rokko_parameters params, const char* key) { 
  return static_cast<rokko::parameters*>(params.ptr)->get<int>(key);
}

int rokko_parameters_get_logicalint(rokko_parameters params, const char* key) { 
  return static_cast<rokko::parameters*>(params.ptr)->get<bool>(key);  // automatically cast bool to int
}

double rokko_parameters_get_double(rokko_parameters params, const char* key) {
  return static_cast<rokko::parameters*>(params.ptr)->get<double>(key);
}

char rokko_parameters_get_char(rokko_parameters params, const char* key) { 
  return static_cast<rokko::parameters*>(params.ptr)->get<char>(key);
}

char* copy_string(std::string const& str) {
  char* p = (char*) malloc( sizeof(char) * (str.size() + 1) );
  str.copy(p, str.size(), 0);
  p[str.size()] = '\0';
  return p;
}


char* rokko_parameters_get_string(struct rokko_parameters params, const char* key) {
  std::string tmp = static_cast<rokko::parameters*>(params.ptr)->get_string(key);
  char* p = copy_string(tmp);
  return p;
}

int rokko_parameters_defined(struct rokko_parameters params, const char* key) {
  return static_cast<int>(static_cast<rokko::parameters*>(params.ptr)->defined(key));
}

int rokko_parameters_get_key_size(struct rokko_parameters params, const char* key) {
  return static_cast<rokko::parameters*>(params.ptr)->get_string(key).size();
}

int rokko_parameters_size(struct rokko_parameters params) {
  return static_cast<rokko::parameters*>(params.ptr)->get_map().size();
}

char** rokko_parameters_keys(struct rokko_parameters params) {
  rokko::parameters* ptr = static_cast<rokko::parameters*>(params.ptr);
  int num_keys = ptr->get_map().size();
  char** keys = (char**) malloc( sizeof(char*) * num_keys );
  int i = 0;
  BOOST_FOREACH(rokko::parameters::value_type const& term, ptr->get_map()) {
    keys[i] = copy_string(term.first);
    ++i;
  }
  return keys;
}

char* rokko_string_i(char** strings, int i) {
  //  std::cout << "i=" << i << " string" << strings[i] << std::endl;
  return strings[i];
}


