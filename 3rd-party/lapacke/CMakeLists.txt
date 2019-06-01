if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h)
  include(FortranCInterface)
  FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h 
                           MACRO_NAMESPACE "LAPACK_"
                           SYMBOL_NAMESPACE "LAPACK_")
  if(NOT FortranCInterface_COMPILED)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/lapacke_mangling.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h COPYONLY)
  endif(NOT FortranCInterface_COMPILED)
endif(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h)

file(GLOB SOURCES src/*.c utils/*.c)

foreach(_func geqpf ggsvd ggsvp)
  try_compile(_try_compile ${CMAKE_CURRENT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/config/check_${_func}.f90
    OUTPUT_VARIABLE LOG)
  message(STATUS "Checking for ${_func} in LAPACK: ${_try_compile}")
  if(NOT _try_compile)
    list(REMOVE_ITEM SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_s${_func}.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_s${_func}_work.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_d${_func}.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_d${_func}_work.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_c${_func}.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_c${_func}_work.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_z${_func}.c
      ${CMAKE_CURRENT_SOURCE_DIR}/src/lapacke_z${_func}_work.c)
  endif(NOT _try_compile)
 endforeach()

add_library(lapacke ${SOURCES})
target_link_libraries(lapacke ${LAPACK_LIBRARIES})
install(TARGETS lapacke ARCHIVE DESTINATION lib LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
install(FILES include/lapacke.h include/lapacke_config.h
  include/lapacke_mangling_with_flags.h include/lapacke_utils.h
  ${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h COMPONENT headers
  DESTINATION include)

if(BUILD_EXAMPLES)
  add_subdirectory(example)
endif(BUILD_EXAMPLES)
