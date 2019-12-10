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

#include <rokko/utility/sort_eigenpairs.hpp>

#include <numeric>
#include <algorithm>

#include <gtest/gtest.h>

template <int MAJOR>
void test(bool ascending) {
  constexpr int num = 10;
  std::vector<int> index(num);
  std::iota(index.begin(), index.end(), 0);
  std::random_shuffle(index.begin(), index.end());
  Eigen::VectorXd eigvals(num), eigvals_sorted(num);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> eigvecs(num, num), eigvecs_sorted(num, num);
  for(int i=0; i<num; ++i) {
    eigvals(i) = 1.0*index[i];
    if (MAJOR == Eigen::RowMajor) {
      eigvecs.row(i).setConstant(eigvals(i));
    } else {
      eigvecs.col(i).setConstant(eigvals(i));
    }
  }
  rokko::sort_eigenpairs(eigvals, eigvecs, eigvals_sorted, eigvecs_sorted, ascending);
  for(int i=0; i<num; ++i){
    std::cout << "dim: " << i << std::endl;
    const double e = ascending ? 1.0 * i : 1.0*(num-i-1);
    if (MAJOR == Eigen::RowMajor) {
      ASSERT_EQ( eigvals_sorted(i), e);
      ASSERT_EQ( eigvecs_sorted(i,0), e);
      ASSERT_EQ( eigvecs_sorted(i,1), e);
    } else {
      ASSERT_EQ( eigvals_sorted(i), e);
      ASSERT_EQ( eigvecs_sorted(0,i), e);
      ASSERT_EQ( eigvecs_sorted(1,i), e);
    }
  }
}

TEST(sort_eigenpairs, matrix_majors) {
  {
    std::cout << "row_major, ascending order\n";
    test<Eigen::RowMajor>(true);
  }
  {
    std::cout << "row_major, descending order\n";
    test<Eigen::RowMajor>(false);
  }
  {
    std::cout << "col_major, ascending order\n";
    test<Eigen::ColMajor>(true);
  }
  {
    std::cout << "col_major, descending order\n";
    test<Eigen::ColMajor>(false);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
