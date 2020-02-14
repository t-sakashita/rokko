/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/factory.hpp>
#include <iostream>
#include <string>
#include <rokko/utility/macro_join.hpp>

class test_base {
public:
  virtual ~test_base() = default;
  virtual void print() = 0;
};

template<typename SOLVER>
class test_wrapper : public test_base {
  using solver_type = SOLVER;

public:
  test_wrapper() = default;
  virtual ~test_wrapper() = default;
  void print() {
    solver_impl_.print();
  }

private:
  solver_type solver_impl_;
};

class test1 {
public:
  test1() = default;
  ~test1() = default;
  void print() {
    std::cout << "print_test1" << std::endl;
  }
};

class test2 {
public:
  test2() = default;
  ~test2() = default;
  void print() {
    std::cout << "print_test2" << std::endl;
  }
};


using test_factory = rokko::factory<test_base>;

template<>
test_factory *test_factory::instance_ = nullptr;

#define ROKKO_REGISTER_TEST(solver, name, priority) \
namespace { namespace ROKKO_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::factory<test_base>; \
  using product = test_wrapper<solver>; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

ROKKO_REGISTER_TEST(test1, "test1", 10)
ROKKO_REGISTER_TEST(test2, "test2", 20)

int main() {
  test_factory::product_pointer_type solver_impl = test_factory::instance()->make_product("test2");
  solver_impl->print();
  std::cerr << "product_names:" << std::endl;
  for(auto name : test_factory::instance()->product_names()) {
    std::cerr << name << std::endl;
  }
  std::cerr << "default_product_name = " << test_factory::instance()->default_product_name() << std::endl;
}
