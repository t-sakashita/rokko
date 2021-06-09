/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>      // for std::istringstream
#include <ctime>

#include <rokko/config.h>

#ifdef ROKKO_HAVE_BOOST
#include <boost/asio.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rokko {

void machine_info(MPI_Comm const& comm = MPI_COMM_WORLD) {
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  std::time_t now = std::time(0);
  std::cout << "num_procs = " << nprocs << std::endl
            #ifdef _OPENMP
            << "num_threads per process = " << omp_get_max_threads() << std::endl
            #endif
#ifdef ROKKO_HAVE_BOOST
            << "hostname = " << boost::asio::ip::host_name() << std::endl
#endif
            << "rokko_version = " << ROKKO_VERSION << std::endl
            << "date = " << ctime(&now);
}


void machine_info(std::ostringstream& oss, MPI_Comm const& comm = MPI_COMM_WORLD) {
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  std::time_t now = std::time(0);
  oss << "num_procs = " << nprocs << std::endl
            #ifdef _OPENMP
	    << "num_threads per process = " << omp_get_max_threads() << std::endl
            #endif
#ifdef ROKKO_HAVE_BOOST
	    << "hostname = " << boost::asio::ip::host_name() << std::endl
#endif
	    << "rokko_version = " << ROKKO_VERSION << std::endl
	    << "date = " << ctime(&now);
}

} // namespace rokko
