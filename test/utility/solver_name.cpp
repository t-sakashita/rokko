/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/utility/solver_name.hpp>

#include <gtest/gtest.h>

TEST(split_solver_name, const_string) {
  std::string library, routine;
  rokko::split_solver_name("lapack:dsyevd", library, routine);
  ASSERT_TRUE(library == "lapack");
  ASSERT_TRUE(routine == "dsyevd");
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
