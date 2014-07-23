include_directories(${PROJECT_SOURCE_DIR}/rokko_c-fortran)
include_directories(${PROJECT_BINARY_DIR}/rokko_c-fortran)

set(SAMPLE frank frank_function_pointer)
foreach(name ${SAMPLE})
  set(target_name sample_c_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(${target_name} rokko_c ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB})
endforeach(name)

