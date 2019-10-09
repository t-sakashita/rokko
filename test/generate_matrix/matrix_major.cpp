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

#include <rokko/eigen3.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gtest/gtest.h>

TEST(matrix, major) {
  int dim = 3;

  // M = 1 2 3
  //     4 5 6
  //     7 8 9

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> M0(dim,dim); // row major
  M0 << 1,2,3,4,5,6,7,8,9;
  std::cout << M0 << std::endl;
  ASSERT_EQ(M0(0,0),1.0);
  ASSERT_EQ(M0(0,1),2.0);
  ASSERT_EQ(M0(1,0),4.0);
  double* ptr0 = &M0(0,0);
  ASSERT_EQ(*(ptr0 + 1), 2.0); // row major
  ASSERT_EQ((M0.row(0))(1), 2.0);
  ASSERT_EQ((M0.col(0))(1), 4.0);

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> M1(dim,dim); // column major
  M1 << 1,2,3,4,5,6,7,8,9;
  std::cout << M1 << std::endl;
  ASSERT_EQ(M1(0,0),1.0);
  ASSERT_EQ(M1(0,1),2.0);
  ASSERT_EQ(M1(1,0),4.0);
  double* ptr1 = &M1(0,0);
  ASSERT_EQ(*(ptr1 + 1), 4.0); // column major
  ASSERT_EQ((M1.row(0))(1), 2.0);
  ASSERT_EQ((M1.col(0))(1), 4.0);

  boost::numeric::ublas::matrix<double> M2(dim,dim); // row major
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M2(i,j) = M0(i,j);
  std::cout << M2 << std::endl;
  ASSERT_EQ(M2(0,0),1.0);
  ASSERT_EQ(M2(0,1),2.0);
  ASSERT_EQ(M2(1,0),4.0);
  double* ptr2 = &M2(0,0);
  ASSERT_EQ(*(ptr2 + 1), 2.0); // row major

  boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> M3(dim,dim); // column major
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M3(i,j) = M0(i,j);
  std::cout << M3 << std::endl;
  ASSERT_EQ(M3(0,0),1.0);
  ASSERT_EQ(M3(0,1),2.0);
  ASSERT_EQ(M3(1,0),4.0);
  double* ptr3 = &M3(0,0);
  ASSERT_EQ(*(ptr3 + 1), 4.0); // column major

  Eigen::MatrixXd M4(dim, dim); // column major
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M4(i,j) = M0(i,j);
  std::cout << M4 << std::endl;
  ASSERT_EQ(M4(0,0),1.0);
  ASSERT_EQ(M4(0,1),2.0);
  ASSERT_EQ(M4(1,0),4.0);
  double* ptr4 = &M4(0,0);
  ASSERT_EQ(*(ptr4 + 1), 4.0); // column major
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
