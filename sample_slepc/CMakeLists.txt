if(BUILD_SLEPC)
  include_directories(${SLEPC_INCLUDE_DIR} ${PETSC_INCLUDE_DIR})
  set(TESTS frank_1proc heisenberg_1proc heisenberg_matrix)

  foreach(name ${TESTS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name}  ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES})
  endforeach(name)
endif(BUILD_SLEPC)
