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

#ifndef ROKKO_JOIN_HPP
#define ROKKO_JOIN_HPP

// based on idea of BOOST_JOIN
#define ROKKO_JOIN(X, Y) ROKKO_DO_JOIN(X, Y)
#define ROKKO_DO_JOIN(X, Y) ROKKO_DO_JOIN2(X,Y)
#define ROKKO_DO_JOIN2(X, Y) X##Y

#endif // ROKKO_JOIN_HPP
