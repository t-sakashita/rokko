/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2010-2013 by Synge Todo <wistaria@comp-phys.org>,
*                            Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>,
*                            Hideyuki Shitara <shitara.hide@jp.fujitsu.com>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTIL_TIMER_HPP
#define ROKKO_UTIL_TIMER_HPP

// ROKKO_ENBALE_TIMER
// ROKKO_ENABLE_TIMER_TRACE
// ROKKO_ENABLE_TIMER_DETAILED

#include <iostream>
#include <string>
#include <boost/config.hpp>

namespace rokko {

class timer_dumb {
public:
  BOOST_STATIC_CONSTANT(int, detailed = (1 << 0));
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
  void summarize(std::ostream& = std::clog) const {}
  void detailed_report(std::size_t = 1, std::ostream& = std::clog) const {}
};

} // end namespace rokko

#ifdef ROKKO_ENABLE_TIMER

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <sys/time.h> /* gettimeofday */
#include <valarray>
#include <vector>
#include <fstream>
#include <map>

#ifdef __linux
#include <sys/types.h>
#include <unistd.h>
#endif

namespace rokko {
namespace detail {

struct clock {
  static double get_time() {
    struct timeval tm;
    gettimeofday(&tm, NULL);
    return tm.tv_sec + tm.tv_usec * 1.0e-6;
  }
};

template<class CLOCK>
class timer_base {
private:
  typedef CLOCK clock_t;
public:
  BOOST_STATIC_CONSTANT(int, detailed = (1 << 0));
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
      boost::throw_exception(std::invalid_argument("empty label"));
    }
    if (id < labels_.size() && !labels_[id].empty()) {
      std::cerr << "Error: duplicated id: " << id << std::endl;
      boost::throw_exception(std::invalid_argument("duplicated id"));
    }
    if (id >= labels_.size()) {
      int old_size = labels_.size();
      int new_size = std::max(2 * labels_.size(), id + 1);
      labels_.resize(new_size);
      std::valarray<double> starts_old = starts_;
      std::valarray<double> counts_old = counts_;
      std::valarray<double> sums_old = sums_;
      starts_.resize(new_size);
      counts_.resize(new_size);
      sums_.resize(new_size);
      for (int i = 0; i < old_size; ++i) {
        starts_[i] = starts_old[i];
        counts_[i] = counts_old[i];
        sums_[i] = sums_old[i];
      }
      #ifdef ROKKO_ENABLE_TIMER_DETAILED
      d_mapping_.resize(new_size, -1);
      #endif
    }
    labels_[id] = label;
    #ifdef ROKKO_ENABLE_TIMER_DETAILED
    if (option & detailed) {
      d_mapping_[id] = d_counts_.size();
      d_counts_.push_back(0);
      d_sums_.push_back(0);
    }
    #endif
    #ifdef ROKKO_ENABLE_TIMER_TRACE
    std::clog << "rokko::timer: registered timer with id = " << id
              << " and label = \"" << label << "\"\n";
    #endif
  }

  void start(std::size_t id, bool fjtool = true) { // hwm: hardware monitor
    #ifdef ROKKO_ENABLE_TIMER_TRACE
    std::clog << "rokko::timer: starting timer with id = " << id << std::endl;
    #endif
    starts_[id] = clock_t::get_time();
  }

  void stop(std::size_t id, bool fjtool = true) {
    double t = clock_t::get_time() - starts_[id];
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

#ifdef __linux
  std::map<std::string, int> vmem_info() const {
    using std::string;

    std::map<string,int> vm_info;
    vm_info.insert(std::pair<string,int>("VmPeak", 0));
    vm_info.insert(std::pair<string,int>("VmHWM", 0));

    // process file name: /proc/${PID}/status
    int pid = static_cast<int>(getppid());
    std::stringstream iss;
    iss << pid;
    std::string filename = string("/proc/")
      + iss.str() + string("/status");

    // open file
    std::ifstream fin(filename.c_str());
    if (fin.fail()) {
      std::cerr << "Can't open the file" << std::endl;
      return std::map<std::string, int>();
    }

    // extract process infomation from procfs
    do {
      string source;
      getline(fin, source);
      for (std::map<string, int>::iterator ip = vm_info.begin(); ip != vm_info.end(); ++ip) {
        if (source.substr(0, (ip->first).size()) == ip->first) {
          string mem_size;
          for (string::size_type i = 0; i != source.size(); ++i) {
            if (isdigit(source[i])) mem_size += source[i];
          }
          ip->second = boost::lexical_cast<int>(mem_size);
        }
      }
    } while (fin.good());
    return vm_info;
  }
#else
  std::map<std::string, int> vmem_info() const {
    std::map<std::string, int> vmem_info;
  }
#endif

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
#ifdef __linux
    std::map<std::string, int> vm = vmem_info();
    for (std::map<std::string, int>::iterator ivm = vm.begin(); ivm != vm.end(); ++ivm) {
      os << "timer: " << ivm->first << "   " << ivm->second << " [kB]" << std::endl;
    }
#endif
    for (int i = 0; i < labels_.size(); ++i) {
      if (counts_[i] > 0) {
        os << boost::format("timer: %5d %-55s %12.3lf %10ld\n")
          % i % labels_[i] % sums_[i] % counts_[i];
      }
    }
  }

  void detailed_report(std::size_t interval = 1, std::ostream& os = std::clog) const {
    #ifdef ROKKO_ENABLE_TIMER_DETAILED
    if (d_count_ == 0) {
      os << "timer: interval = " << interval << std::endl;
    }
    if ((d_count_ + 1) % interval == 0) {
      for (int i = 0; i < labels_.size(); ++i) {
        if (d_mapping_[i] >= 0) {
          os << boost::format("detail: %d %d %.3f %d\n")
            % d_count_ % i % d_sums_[d_mapping_[i]] % d_counts_[d_mapping_[i]];
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
  std::valarray<double> starts_, counts_, sums_;
  #ifdef ROKKO_ENABLE_TIMER_DETAILED
  mutable int d_count_;
  std::vector<int> d_mapping_;
  mutable std::vector<int> d_counts_;
  mutable std::vector<double> d_sums_;
  #endif
  std::map<std::string, int> vm_info_;
};

} // namespace detail

typedef detail::timer_base<detail::clock> timer;

} // namespace rokko

#else

namespace rokko {

typedef timer_dumb timer;

} // namespace rokko

#endif

#endif // ROKKO_UTIL_TIMER_HPP
