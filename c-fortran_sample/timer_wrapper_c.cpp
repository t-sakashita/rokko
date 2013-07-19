/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*    
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/limits.hpp>
#include <rokko/utility/timer.hpp>
#include "timer_wrapper_c.h"

namespace rokko {
  void* initialize_timer(){ 
    
    return static_cast<void*>(new timer_dumb());
  }
  void delete_timer(void* timer){ 
    timer_dumb* timer_ =static_cast<timer_dumb*>(timer);
    delete timer_;
  }
  void timer_start(void* timer, unsigned int id){ 
    timer_dumb* timer_ = static_cast<timer_dumb*>(timer);
    timer_->start(id);
  }
  void timer_stop(void* timer, unsigned int id){ 
    timer_dumb* timer_ = static_cast<timer_dumb*>(timer);
    timer_->stop(id);
  }
  void timer_registrate(void* timer, unsigned int id, char * label){  
    timer_dumb* timer_ = static_cast<timer_dumb*>(timer);
    std::string label_(label);
    timer_->registrate(id, label_);
  }
  double timer_get_count(void* timer, unsigned int id){
    timer_dumb* timer_ = static_cast<timer_dumb*>(timer);
    return timer_->get_count(id);
  }    
  int timer_get_average(void* timer, unsigned int id){
    timer_dumb* timer_ = static_cast<timer_dumb*>(timer);
    return timer_->get_average(id);
  }    



}
