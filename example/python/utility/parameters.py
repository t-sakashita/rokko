#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from rokko import *

params = rokko_parameters()
params.set("PE", 2.3)

b = int(2)
params.set("ABCD", b)

c = "pppppp"
params.set("STR", c)

b = params.get("ABCD")
c = params.get("STR")

print(params.get("ABCD"))
print(params.get_string("ABCD"))
print "b=", b

dic = params.dict()
print(dic["STR"])

print(params.keys())

#params.clear()
