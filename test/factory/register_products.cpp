/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/factory.hpp>


class test_base {
public:
  virtual ~test_base() {};
  virtual void print() = 0;
};

template<typename SOLVER>
class test_wrapper : public test_base {
public:
  test_wrapper() {}
  virtual ~test_wrapper() {}
  void print() {}
};

class test1 {
public:
  test1() {}
  ~test1() {}
  void print() {}
};

class test2 {
public:
  test2() {}
  ~test2() {}
  void print() {}
};


typedef rokko::factory<test_base> test_solver_factory;

template<>
test_solver_factory *test_solver_factory::instance_ = 0;

#define ROKKO_REGISTER_TEST(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  typedef rokko::factory<test_base> factory; \
  typedef test_wrapper<solver> product; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

ROKKO_REGISTER_TEST(test1, "test1", 10)
ROKKO_REGISTER_TEST(test2, "test2", 20)

int main()
{
}
