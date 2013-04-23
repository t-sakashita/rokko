#ifndef ROKKO_CORE_ELEMENTAL_H
#define ROKKO_CORE_ELEMENTAL_H

#include <rokko/elemental/elemental.hpp>

#include <boost/noncopyable.hpp>

namespace rokko {

template<typename T>
void Initialize(int& argc, char**& argv)
{
}

template<typename T>
void Finalize()
{
}

template<>
void Initialize<elemental>(int& argc, char**& argv)
{
  elem::Initialize(argc, argv);
}

template<>
void Finalize<elemental>()
{
  elem::Finalize();
}

} // namespace rokko

#endif // ROKKO_CORE_ELEMENTAL_H
