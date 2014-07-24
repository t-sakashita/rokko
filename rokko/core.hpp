#ifndef ROKKO_CORE_HPP
#define ROKKO_CORE_HPP

#include <boost/noncopyable.hpp>

namespace rokko {

template<typename T>
void initialize(int& argc, char**& argv);

template<typename T>
void finalize();


} // namespace rokko

#endif // ROKKO_CORE_HPP

