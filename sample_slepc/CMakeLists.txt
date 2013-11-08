if(BUILD_SLEPC)
  include_directories(${SLEPC_INCLUDE_DIR} ${PETSC_INCLUDE_DIR})
  set(TESTS frank_sparse_1proc heisenberg_matfree_1proc heisenberg_matfree_mpi heisenberg_sparse_mpi xyz_matfree_mpi)

  foreach(name ${TESTS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name}  ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES})
  endforeach(name)
endif(BUILD_SLEPC)
