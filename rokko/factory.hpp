/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/noncopyable.hpp>
#include <rokko/utility/macro_join.hpp>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <map>
#include <string>
#include <vector>

namespace rokko {

template<typename BASE, class... Types>
class factory : private rokko::noncopyable<factory<BASE>> {
private:
  class abstract_creator {
  public:
    virtual ~abstract_creator() = default;
    virtual std::shared_ptr<BASE> create(Types... args) const = 0;
  };
  template<typename PRODUCT>
  class creator : public abstract_creator {
  public:
    virtual ~creator() = default;
    std::shared_ptr<BASE> create(Types... args) const override {
      return std::static_pointer_cast<BASE>(std::make_shared<PRODUCT>(args...));
    }
  };

public:
  using product_pointer_type = std::shared_ptr<BASE>;

private:
  using creator_pointer_type = std::shared_ptr<abstract_creator>;
  using creator_map_type = std::map<std::string, creator_pointer_type>;
public:
  factory() : largest_priority_(0) {}
  static product_pointer_type make_product(std::string const& name = "", Types... args) {
    factory* f = factory::instance();
    if (name == "") {
      if (f->default_product_ != "") {
        return f->make_creator(f->default_product_)->create(args...);
      } else {
        std::cerr << "Error: default product is not defined\n";
        throw std::runtime_error("factory::make_product()");
      }
    } else {
      return f->make_creator(name)->create(args...);
    }
  }
  template<typename PRODUCT>
  bool register_creator(std::string const& name, int priority = 0) {
    bool isnew = (creators_.find(name) == creators_.cend());
    creators_[name] = std::static_pointer_cast<abstract_creator>(std::make_shared<creator<PRODUCT>>());
    if (priority >= largest_priority_) {
      largest_priority_ = priority;
      default_product_ = name;
    }
    return isnew;
  }
  bool unregister_creator(std::string const& name); // to be implemented
  static std::vector<std::string> product_names() {
    factory* f = factory::instance();
    std::vector<std::string> retvec;
    for (typename creator_map_type::const_iterator it = f->creators_.cbegin();
         it != f->creators_.cend(); ++it) {
      retvec.emplace_back(it->first);
    }
    return retvec;
  }

  static const std::string& default_product_name() {
    return instance()->default_product_;
  }
  static factory* instance() {
    if (!instance_) instance_ = new factory;
    return instance_;
  }
protected:
  creator_pointer_type make_creator(std::string const& name) const {
    typename creator_map_type::const_iterator itr = creators_.find(name);
    if (itr == creators_.cend() || itr->second == nullptr) {
      std::cerr << "Error: unknown product: \"" << name << "\" (registered products: ";
      for (typename creator_map_type::const_iterator itr = creators_.cbegin();
           itr != creators_.cend(); ++itr) {
        if (itr != creators_.cbegin()) std::cerr << ", ";
        std::cerr << "\"" << itr->first << "\"";
      }
      std::cerr << ")\n";
      throw std::runtime_error("factory::make_creator()");
    }
    return itr->second;
  }
      
private:
  static inline factory* instance_;
  creator_map_type creators_;
  int largest_priority_;
  std::string default_product_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PRODUCT(base, product, name, priority) \
  namespace { namespace ROKKO_JOIN(product_register, __LINE__) { struct register_caller { register_caller() { rokko::factory::instance()->register_creator<product>(name, priority); } } caller; } }
