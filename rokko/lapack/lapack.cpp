#ifndef ROKKO_LAPACK_CPP
#define ROKKO_LAPACK_CPP

#include <rokko/serial_solver_factory.hpp>
#include <rokko/lapack/core.hpp>

ROKKO_REGISTER_SERIAL_SOLVER(rokko::lapack::solver<rokko::lapack::dsyev>, "lapack")

#endif
