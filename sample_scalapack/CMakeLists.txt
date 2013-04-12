cmake_minimum_required(VERSION 2.8)

project(frank NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# enable C and C++ compilers
enable_language(CXX)

find_package(MPI)

# find Eigen3 Library
find_package(Eigen3)

include_directories(BEFORE ${EIGEN3_INCLUDE_DIR})

set(TESTS  diagonalize frank_old product diagonalize2)

foreach(name ${TESTS})
  add_executable(${name}  ${name}.cpp)
  target_link_libraries(${name}  "-lmkl_scalapack_lp64" "-lmkl_blacs_intelmpi_lp64" "-mkl=sequential")
  target_link_libraries(${name} ${MPI_CXX_LIBRARIES})
  #add_alps_test(${name})
endforeach(name)



