/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/localized_matrix.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define BOOST_TEST_MODULE test_matrix_major
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_matrix_major) {
  int dim = 3;

  // M = 1 2 3
  //     4 5 6
  //     7 8 9

  rokko::localized_matrix<double, rokko::matrix_row_major> M0(dim,dim); // row major
  M0 << 1,2,3,4,5,6,7,8,9;
  std::cout << M0 << std::endl;
  BOOST_CHECK_EQUAL(M0(0,0),1.0);
  BOOST_CHECK_EQUAL(M0(0,1),2.0);
  BOOST_CHECK_EQUAL(M0(1,0),4.0);
  double* ptr0 = &M0(0,0);
  BOOST_CHECK_EQUAL(*(ptr0 + 1), 2.0); // row major
  BOOST_CHECK_EQUAL((M0.row(0))(1), 2.0);
  BOOST_CHECK_EQUAL((M0.col(0))(1), 4.0);

  rokko::localized_matrix<double, rokko::matrix_col_major> M1(dim,dim); // column major
  M1 << 1,2,3,4,5,6,7,8,9;
  std::cout << M1 << std::endl;
  BOOST_CHECK_EQUAL(M1(0,0),1.0);
  BOOST_CHECK_EQUAL(M1(0,1),2.0);
  BOOST_CHECK_EQUAL(M1(1,0),4.0);
  double* ptr1 = &M1(0,0);
  BOOST_CHECK_EQUAL(*(ptr1 + 1), 4.0); // column major
  BOOST_CHECK_EQUAL((M1.row(0))(1), 2.0);
  BOOST_CHECK_EQUAL((M1.col(0))(1), 4.0);

  boost::numeric::ublas::matrix<double> M2(dim,dim); // row major
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M2(i,j) = M0(i,j);
  std::cout << M2 << std::endl;
  BOOST_CHECK_EQUAL(M2(0,0),1.0);
  BOOST_CHECK_EQUAL(M2(0,1),2.0);
  BOOST_CHECK_EQUAL(M2(1,0),4.0);
  double* ptr2 = &M2(0,0);
  BOOST_CHECK_EQUAL(*(ptr2 + 1), 2.0); // row major

  boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> M3(dim,dim); // column major
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M3(i,j) = M0(i,j);
  std::cout << M3 << std::endl;
  BOOST_CHECK_EQUAL(M3(0,0),1.0);
  BOOST_CHECK_EQUAL(M3(0,1),2.0);
  BOOST_CHECK_EQUAL(M3(1,0),4.0);
  double* ptr3 = &M3(0,0);
  BOOST_CHECK_EQUAL(*(ptr3 + 1), 4.0); // column major

  Eigen::MatrixXd M4(dim, dim); // column major
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M4(i,j) = M0(i,j);
  std::cout << M4 << std::endl;
  BOOST_CHECK_EQUAL(M4(0,0),1.0);
  BOOST_CHECK_EQUAL(M4(0,1),2.0);
  BOOST_CHECK_EQUAL(M4(1,0),4.0);
  double* ptr4 = &M4(0,0);
  BOOST_CHECK_EQUAL(*(ptr4 + 1), 4.0); // column major
}
