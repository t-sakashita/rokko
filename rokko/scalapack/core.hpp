#ifndef ROKKO_SCALAPACK_CORE_HPP
#define ROKKO_SCALAPACK_CORE_HPP

#include <rokko/scalapack/scalapack.hpp>

#include <boost/noncopyable.hpp>

#include <rokko/core.hpp>

namespace rokko {

template<>
void initialize<rokko::scalapack>(int& argc, char**& argv)
{
}

template<>
void finalize<rokko::scalapack>()
{
}

} // namespace rokko

#endif // ROKKO_SCALAPACK_CORE_HPP

