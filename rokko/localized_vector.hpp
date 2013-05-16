#ifndef ROKKO_LOCALIZED_VECTOR_HPP
#define ROKKO_LOCALIZED_VECTOR_HPP

#include <Eigen/Dense>

namespace rokko {

struct localized_vector : public Eigen::VectorXd {
public:
  localized_vector() : Eigen::VectorXd() {};
  localized_vector(int size) : Eigen::VectorXd(size) {};
};


} // namespace rokko

#endif // ROKKO_LOCALIZED_VECTOR_HPP
