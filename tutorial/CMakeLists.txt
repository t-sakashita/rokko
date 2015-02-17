#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2014 by Synge Todo <wistaria@comp-phys.org>
# Copyright (C) 2015 by Yuichi Motoyama <y-motoyama@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

install(DIRECTORY . DESTINATION tutorial COMPONENT tutorial
  FILES_MATCHING PATTERN "*.txt"
                 PATTERN "*.md"
                 PATTERN "*.C"
                 PATTERN "*.cpp"
                 PATTERN "*.hpp"
                 PATTERN "*.c"
                 PATTERN "*.h"
                 PATTERN "*.f"
                 PATTERN "*.f90"
                 PATTERN "*.F90"
                 PATTERN "*.ip"
                 PATTERN "*.op"
                 )
