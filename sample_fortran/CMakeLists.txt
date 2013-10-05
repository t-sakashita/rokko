include_directories(${PROJECT_SOURCE_DIR}/rokko_c-fortran)
include_directories(${PROJECT_BINARY_DIR}/rokko_c-fortran)

set(SAMPLE frank frank_function frank_array frank_array_all_gather set_distributed_matrix)
#frank_function_array
foreach(name ${SAMPLE})
  set(target_name sample_fortran_${name})
  add_executable(${target_name} ${name}.f90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko_fortran ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB})
  SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach(name)

