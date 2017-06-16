add_subdirectory(cblas)
add_subdirectory(lapacke)

if(BUILD_CXX_CORE)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/3rd-party/eigen3/cmake)
  set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
  add_subdirectory(eigen3/Eigen)
endif(BUILD_CXX_CORE)
