#ifndef ROKKO_EIGEN_S_CORE_HPP
#define ROKKO_EIGEN_S_CORE_HPP

#include <rokko/eigen_s/eigen_s.hpp>

#include <boost/noncopyable.hpp>

#include <rokko/core.hpp>

namespace rokko {

template<>
void Initialize<rokko::eigen_s>(int& argc, char**& argv)
{
}

template<>
void Finalize<rokko::eigen_s>()
{
}

} // namespace rokko

#endif // ROKKO_EIGEN_S_CORE_HPP

