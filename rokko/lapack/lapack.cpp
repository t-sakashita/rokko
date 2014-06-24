#ifndef ROKKO_LAPACK_CPP
#define ROKKO_LAPACK_CPP

#include <rokko/solver_factory.hpp>
#include <rokko/lapack/core.hpp>

ROKKO_REGISTER_SERIAL_DENSE_SOLVER(rokko::lapack::solver<rokko::lapack::dsyev>, "lapack", 20)

#endif
