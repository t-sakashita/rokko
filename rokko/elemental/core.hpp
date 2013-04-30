#ifndef ROKKO_CORE_ELEMENTAL_H
#define ROKKO_CORE_ELEMENTAL_H

#include <rokko/elemental/elemental.hpp>

#include <boost/noncopyable.hpp>

namespace rokko {

template<typename T>
void initialize(int& argc, char**& argv)
{
}

template<typename T>
void finalize()
{
}

template<>
void initialize<elemental>(int& argc, char**& argv)
{
  elem::Initialize(argc, argv);
}

template<>
void finalize<elemental>()
{
  elem::Finalize();
}

} // namespace rokko

#endif // ROKKO_CORE_ELEMENTAL_H
