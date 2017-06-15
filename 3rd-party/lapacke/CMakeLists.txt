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
add_library(lapacke ${SOURCES})
if (NOT STATIC_SCALAPACK)
  target_link_libraries(lapacke ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})
endif (NOT STATIC_SCALAPACK)
install(TARGETS lapacke ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
install(FILES include/lapacke.h include/lapacke_config.h include/lapacke_mangling_with_flags.h include/lapacke_utils.h ${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h COMPONENT headers DESTINATION include)

add_subdirectory(example)
