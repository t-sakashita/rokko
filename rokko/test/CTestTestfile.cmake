# CMake generated Testfile for 
# Source directory: /home/sakashita/development/rokko/test
# Build directory: /home/sakashita/development/rokko/rokko/test
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(scatter "/home/issp/usr/bin/cmake" "-Dcmd=scatter" "-Dsourcedir=/home/sakashita/development/rokko/test" "-Dbinarydir=/home/sakashita/development/rokko/rokko/test" "-Ddllexedir=/home/sakashita/development/rokko/rokko/bin" "-Dinput=scatter" "-Doutput=scatter" "-P" "/home/sakashita/development/rokko/config/run_test.cmake")
ADD_TEST(gather "/home/issp/usr/bin/cmake" "-Dcmd=gather" "-Dsourcedir=/home/sakashita/development/rokko/test" "-Dbinarydir=/home/sakashita/development/rokko/rokko/test" "-Ddllexedir=/home/sakashita/development/rokko/rokko/bin" "-Dinput=gather" "-Doutput=gather" "-P" "/home/sakashita/development/rokko/config/run_test.cmake")
