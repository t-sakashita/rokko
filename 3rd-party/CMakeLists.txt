if(BUILD_CBLAS)
  add_subdirectory(cblas)
endif(BUILD_CBLAS)

if(BUILD_LAPACKE)
  add_subdirectory(lapacke)
endif(BUILD_LAPACKE)

if(BUILD_BASE_CXX)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/3rd-party/eigen3/cmake)
  set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
  add_subdirectory(eigen3/Eigen)
endif(BUILD_BASE_CXX)
