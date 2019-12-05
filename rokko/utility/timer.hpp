/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2010-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTIL_TIMER_HPP
#define ROKKO_UTIL_TIMER_HPP

// You may define the following macros:
//   #define ROKKO_DISABLE_TIMER
//   #define ROKKO_ENABLE_TIMER_TRACE
//   #define ROKKO_ENABLE_TIMER_DETAILED

#include <rokko/utility/string_format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/noncopyable.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <chrono>
#include <valarray>
#include <vector>

#ifdef __linux
#include <sys/types.h>
#include <unistd.h>
#endif

namespace rokko {

struct timer_id {
  constexpr static int solver_construct = 1;
  constexpr static int solver_initialize = 2;
  constexpr static int solver_finalize = 3;
  constexpr static int diagonalize_initialize = 4;
  constexpr static int diagonalize_diagonalize = 5;
  constexpr static int diagonalize_finalize = 6;
};

namespace detail {

class timer_base {
public:
  constexpr static int detailed = 1 << 0;
  timer_base() {
    #ifdef ROKKO_ENABLE_TIMER_DETAILED
    d_count_ = 0;
    #endif
  }
  void clear() {
    counts_ = 0;
    sums_ = 0;
    #ifdef ROKKO_TRACE_TIMER
    std::clog << "rokko::timer: cleared\n";
    #endif
  }
  void registrate(std::size_t id, std::string const& label, int option = 0) {
    if (label.empty()) {
      std::cerr << "Error: empty label\n";
      throw std::invalid_argument("empty label");
    }
    if (id < labels_.size() && !labels_[id].empty()) {
      std::cerr << "Error: duplicated id: " << id << std::endl;
      throw std::invalid_argument("duplicated id");
    }
    if (id >= labels_.size()) {
      std::size_t old_size = labels_.size();
      std::size_t new_size = std::max(2 * labels_.size(), id + 1);
      labels_.resize(new_size);
      std::valarray<std::chrono::system_clock::time_point> starts_old = starts_;
      std::valarray<double> counts_old = counts_;
      std::valarray<double> sums_old = sums_;
      starts_.resize(new_size);
      counts_.resize(new_size);
      sums_.resize(new_size);
      for (std::size_t i = 0; i < old_size; ++i) {
        starts_[i] = starts_old[i];
        counts_[i] = counts_old[i];
        sums_[i] = sums_old[i];
      }
      #ifdef ROKKO_ENABLE_TIMER_DETAILED
      d_mapping_.resize(new_size, -1);
      #endif
    }
    labels_[id] = label;
    if (option & detailed) {
      #ifdef ROKKO_ENABLE_TIMER_DETAILED
      d_mapping_[id] = d_counts_.size();
      d_counts_.push_back(0);
      d_sums_.push_back(0);
      #endif
    }
    #ifdef ROKKO_ENABLE_TIMER_TRACE
    std::clog << "rokko::timer: registered timer with id = " << id
              << " and label = \"" << label << "\"\n";
    #endif
  }

  void start(std::size_t id) { // hwm: hardware monitor
    #ifdef ROKKO_ENABLE_TIMER_TRACE
    std::clog << "rokko::timer: starting timer with id = " << id << std::endl;
    #endif
    starts_[id] = std::chrono::system_clock::now();
  }

  void stop(std::size_t id) {
    std::chrono::system_clock::time_point end_ = std::chrono::system_clock::now();
    double t = 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(end_ - starts_[id]).count();
    counts_[id] += 1;
    sums_[id] += t;
    #ifdef ROKKO_ENABLE_TIMER_DETAILED
    if (d_mapping_[id] >= 0) {
      d_counts_[d_mapping_[id]] += 1;
      d_sums_[d_mapping_[id]] += t;
    }
    #endif
    #ifdef ROKKO_ENABLE_TIMER_TRACE
    std::clog << "rokko::timer: stopping timer with id = " << id << std::endl;
    #endif
  }
  std::vector<std::string> const& get_labels() const { return labels_; }
  std::valarray<double> const& get_counts() const { return counts_; }
  std::valarray<double> const& get_measurements() const { return sums_; }

  bool has(std::size_t id) { return (id < labels_.size() && !labels_[id].empty()); }
  std::string const& get_label(std::size_t id) const { return labels_[id]; }
  double get_count(std::size_t id) const { return counts_[id]; }
  double get_measurement(std::size_t id) const { return sums_[id]; }
  double get_average(std::size_t id) const {
    return (counts_[id] > 0 ? sums_[id] / counts_[id] : 0.0);
  }

  std::map<std::string, int> vmem_info() const {
    std::map<std::string, int> vm_info;
// #ifdef __linux
//     vm_info.insert(std::pair<std::string, int>("VmPeak", 0));
//     vm_info.insert(std::pair<std::string, int>("VmHWM", 0));

//     // process file name: /proc/${PID}/status
//     int pid = static_cast<int>(getpid());
//     std::string filename = std::string("/proc/") + boost::lexical_cast<std::string>(pid) +
//       "/status";

//     // open file
//     std::ifstream fin(filename);
//     if (fin.fail()) {
//       std::cerr << "Can't open the file" << std::endl;
//       return vm_info;
//     }

//     // extract process infomation from procfs
//     do {
//       std::string source;
//       getline(fin, source);
//       for (std::map<std::string, int>::iterator ip = vm_info.begin(); ip != vm_info.end(); ++ip) {
//         if (source.substr(0, (ip->first).size()) == ip->first) {
//           std::string mem_size;
//           for (std::string::size_type i = 0; i != source.size(); ++i) {
//             if (isdigit(source[i])) mem_size += source[i];
//           }
//           ip->second = boost::lexical_cast<int>(mem_size);
//         }
//       }
//     } while (fin.good());
// #endif
    return vm_info;
  }

  void summarize(std::ostream& os = std::clog) const {
    os << "timer: enabled\n"
#ifdef ROKKO_ENABLE_TIMER_TRACE
       << "timer: trace = enabled\n"
#else
       << "timer: trace = disabled\n"
#endif
#ifdef ROKKO_ENABLE_TIMER_DETAILED
       << "timer: detailed report = enabled"
#else
       << "timer: detailed report = disabled"
#endif
       << std::endl;
    // std::map<std::string, int> vm = vmem_info();
    // for (std::map<std::string, int>::iterator ivm = vm.begin(); ivm != vm.end(); ++ivm) {
    //   os << "timer: " << ivm->first << "   " << ivm->second << " [kB]" << std::endl;
    // }
    for (std::size_t i = 0; i < labels_.size(); ++i) {
      if (counts_[i] > 0) {
        os << rokko::format("timer: %5d %-55s %12.3lf %10ld\n",
                            i, labels_[i].c_str(), sums_[i], counts_[i]);
      }
    }
  }

  void detailed_report(std::size_t interval = 1, std::ostream& os = std::clog) const {
    #ifdef ROKKO_ENABLE_TIMER_DETAILED
    if (d_count_ == 0) {
      os << "timer: interval = " << interval << std::endl;
    }
    if ((d_count_ + 1) % interval == 0) {
      for (std::size_t i = 0; i < labels_.size(); ++i) {
        if (d_mapping_[i] >= 0) {
          os << boost::format("detail: %d %d %.3f %d\n",
                              d_count_, i, d_sums_[d_mapping_[i]], d_counts_[d_mapping_[i]]);
        }
      }
    }
    ++d_count_;
    std::fill(d_counts_.begin(), d_counts_.end(), 0);
    std::fill(d_sums_.begin(), d_sums_.end(), 0);
    os << std::flush;
    #endif
  }

protected:
  std::vector<std::string> labels_;
  std::valarray<std::chrono::system_clock::time_point> starts_;
  std::valarray<double> counts_, sums_;
  #ifdef ROKKO_ENABLE_TIMER_DETAILED
  mutable int d_count_;
  std::vector<int> d_mapping_;
  mutable std::vector<int> d_counts_;
  mutable std::vector<double> d_sums_;
  #endif
  std::map<std::string, int> vm_info_;
};

class timer_dumb {
public:
  constexpr static int detailed = 0;
  timer_dumb() {}
  void clear() {}
  void registrate(std::size_t, std::string const&, int = 0) {}
  void start(std::size_t) {}
  void stop(std::size_t) {}
  bool has(std::size_t) { return false; }
  std::string get_label(std::size_t) const { return std::string(); }
  double get_count(std::size_t) const { return 0; }
  double get_measurement(std::size_t) const { return 0; }
  double get_average(std::size_t) const { return 0; }
  std::map<std::string, int> vmem_info() const { return std::map<std::string, int>(); }
  void summarize(std::ostream& = std::clog) const {}
  void detailed_report(std::size_t = 1, std::ostream& = std::clog) const {}
};

} // namespace detail

#ifndef ROKKO_DISABLE_TIMER
using timer = detail::timer_base;
#else
using timer = detail::timer_dumb;
#endif

class global_timer : private boost::noncopyable {
public:
  constexpr static int detailed = rokko::timer::detailed;
  static void clear() { instance()->clear(); }
  static void registrate(std::size_t id, std::string const& label, int option = 0) {
    instance()->registrate(id, label, option);
  }
  static void start(std::size_t id) { instance()->start(id); }
  static void stop(std::size_t id) { instance()->stop(id); }
  static bool has(std::size_t id) { return instance()->has(id); }
  static std::string get_label(std::size_t id) { return instance()->get_label(id); }
  static double get_count(std::size_t id) { return instance()->get_count(id); }
  static double get_measurement(std::size_t id) { return instance()->get_measurement(id); }
  static double get_average(std::size_t id) { return instance()->get_average(id); }
  static std::map<std::string, int> vmem_info() { return instance()->vmem_info(); }
  static void summarize(std::ostream& os = std::clog) { instance()->summarize(os); }
  static rokko::timer* instance() {
    if (!instance_) instance_ = new rokko::timer;
    return instance_;
  }
private:
  static rokko::timer* instance_;
};

} // namespace rokko

#endif // ROKKO_UTIL_TIMER_HPP
