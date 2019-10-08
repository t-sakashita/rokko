/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_VECTOR_HPP
#define ROKKO_DISTRIBUTED_VECTOR_HPP

namespace rokko {

template<typename T>
class distributed_vector {
public:
  using value_type = T;
  distributed_vector() {}
  distributed_vector(int n_global, int begin_local, int end_local) {
    initialize(n_global, begin_local, end_local);
  }
  void initialize(int n_global, int begin_local, int end_local) {
    n_global_ = n_global;
    offset_ = begin_local;
    if (storage_.size() < end_local - begin_local) storage_.resize(end_local - begin_local);
  }

  value_type* get_storage() { return &storage_[0]; }
  const value_type* get_storage() const { return &storage_[0]; }
  int size() const { return n_global_; }
  int size_local() const { return storage_.size(); }

  bool is_gindex(int gi) const {
    return (gi >= offset_ && gi < offset_ + storage_.size());
  }

  void set_local(int li, value_type value) { storage_[li] = value; }
  void update_local(int li, value_type value) { storage_[li] += value; }
  value_type get_local(int li) const { return storage_[li]; }

  void set_global(int gi, value_type value) {
    if (is_gindex(gi)) set_local(gi - offset_, value);
  }
  void update_global(int gi, value_type value) {
    if (is_gindex(gi)) update_local(gi - offset_, value);
  }
  value_type get_global(int gi) const { return get_local(gi - offset_); }
  value_type get_global_checked(int gi) const {
    if (is_gindex(gi)) {
      return get_global(gi);
    } else {
      throw std::out_of_range("element not on this process.");
    }
  }

  void set_zeros() { std::fill(storage_.begin(), storage_.end(), 0); }

private:
  int n_global_, offset_;
  std::vector<value_type> storage_;
};

} // namespace rokko

#endif // ROKKO_DISTRIBUTED_VECTOR_HPP
