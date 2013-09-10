# program using SLEPc
message(STATUS "Display PETSC_INCLUDE_DIRS=${PETSC_INCLUDE_DIR}.")
message(STATUS "Display PETSC_LIBRARIES=${PETSC_LIBRARIES}.")
message(STATUS "Display SLEPC_INCLUDE_DIRS=${SLEPC_INCLUDE_DIR}.")
message(STATUS "Display SLEPC_LIBRARIES=${SLEPC_LIBRARIES}.")

include_directories(${SLEPC_INCLUDE_DIR} ${PETSC_INCLUDE_DIR})

set(TESTS frank_1proc heisenberg_1proc heisenberg_matrix)
# heisenberg_parallel)

foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name}  ${SLEPC_LIBRARIES} ${PETSC_LIBRARIES})
endforeach(name)
