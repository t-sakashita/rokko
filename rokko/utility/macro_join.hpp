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

#pragma once

// based on idea of BOOST_JOIN
#define ROKKO_JOIN(X, Y) ROKKO_DO_JOIN(X, Y)
#define ROKKO_DO_JOIN(X, Y) ROKKO_DO_JOIN2(X,Y)
#define ROKKO_DO_JOIN2(X, Y) X##Y
