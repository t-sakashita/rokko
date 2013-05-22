add_subdirectory(lapacke)
add_subdirectory(eigen_s)
add_subdirectory(eigen_sx)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/3rd-party/eigen3/cmake)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
add_subdirectory(eigen3/Eigen)
