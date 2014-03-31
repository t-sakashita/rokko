/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Yuichi Motoyama <y-motoyama@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/utility/sort_eigenpairs.hpp>

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

#define BOOST_TEST_MODULE test_sort_eigenpairs
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#define make_test(major, ascending) \
  int num = 10;\
  std::vector<int> index;\
  std::copy(boost::counting_iterator<int>(0),\
            boost::counting_iterator<int>(num),\
            back_inserter(index));\
  std::random_shuffle(index.begin(), index.end());\
  rokko::localized_vector  eigvals(num), eigvals_sorted(num);\
  rokko::localized_matrix< major > eigvecs(num, num), eigvecs_sorted(num, num);\
  for(int i=0; i<num; ++i){\
    eigvals(i) = 1.0*index[i];\
    for(int j=0; j<num; ++j){\
      if(eigvecs.is_row_major()){\
        eigvecs(i,j) = eigvals(i);\
      }else{\
        eigvecs(j,i) = eigvals(i);\
      }\
    }\
  }\
  rokko::sort_eigenpairs(eigvals, eigvecs, eigvals_sorted, eigvecs_sorted, ascending);\
  for(int i=0; i<num; ++i){\
    std::cout << "dim: " << i << std::endl;\
    double e = ascending ? 1.0 * i : 1.0*(num-i-1);\
    if(eigvecs.is_row_major()){\
      BOOST_CHECK_EQUAL( eigvals_sorted(i), e);\
      BOOST_CHECK_EQUAL( eigvecs_sorted(i,0), e);\
      BOOST_CHECK_EQUAL( eigvecs_sorted(i,1), e);\
    }else{\
      BOOST_CHECK_EQUAL( eigvals_sorted(i), e);\
      BOOST_CHECK_EQUAL( eigvecs_sorted(0,i), e);\
      BOOST_CHECK_EQUAL( eigvecs_sorted(1,i), e);\
    }\
  }

BOOST_AUTO_TEST_CASE(test_sort_eigenpairs) {
  {
    std::cout << "row_major, ascending order\n";
    make_test(rokko::matrix_row_major, true)
  }
  {
    std::cout << "row_major, descending order\n";
    make_test(rokko::matrix_row_major, false)
  }
  {
    std::cout << "col_major, ascending order\n";
    make_test(rokko::matrix_col_major, true)
  }
  {
    std::cout << "col_major, descending order\n";
    make_test(rokko::matrix_col_major, false)
  }
}

#undef make_test
