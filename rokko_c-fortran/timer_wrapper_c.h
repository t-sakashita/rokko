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

#ifdef __cplusplus
namespace rokko {
extern "C"{
#endif

  void* initialize_timer();
  void delete_timer(void*);
  void timer_start(void*, unsigned int);
  void timer_stop(void*, unsigned int);
  void timer_registrate(void*, unsigned int, char*);
  double timer_get_count(void*, unsigned int);
  int timer_get_average(void*, unsigned int);

#ifdef __cplusplus
}
}
#endif
