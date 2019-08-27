#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pyrokko import *

params = parameters()
params.set("PE", 2.3)
params.set("ABCD", 2)
params.set("STR", "pppppp")
params.set("flag", True)

b = params.get("ABCD")
c = params.get("STR")

print("ABCD: ", params.get("ABCD"))
print("flag: ", params.get("flag"))
print("b=", b)

dic = params.dict
print(dic["STR"])

print(params.keys)

params.clear()
