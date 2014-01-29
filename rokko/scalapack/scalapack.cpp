#ifndef ROKKO_SCALAPACK_CPP
#define ROKKO_SCALAPACK_CPP

#include <rokko/solver_factory.hpp>
#include <rokko/scalapack/core.hpp>

ROKKO_REGISTER_SOLVER(rokko::scalapack::solver<rokko::scalapack::pdsyev>, "scalapack")

#endif
