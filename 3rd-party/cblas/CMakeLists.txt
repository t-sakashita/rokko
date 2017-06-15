if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/include/cblas_mangling.h)
  include(FortranCInterface)
  FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/cblas_mangling.h 
                           MACRO_NAMESPACE "BLAS_"
                           SYMBOL_NAMESPACE "BLAS_")
  if(NOT FortranCInterface_COMPILED)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/cblas_mangling.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/cblas_mangling.h COPYONLY)
  endif(NOT FortranCInterface_COMPILED)
endif(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/include/cblas_mangling.h)

file(GLOB SOURCES src/*.c src/*.f)
add_library(cblas ${SOURCES})
target_link_libraries(cblas ${BLAS_LIBRARIES})
install(TARGETS cblas ARCHIVE DESTINATION lib LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(FILES include/cblas.h include/cblas_f77.h
        ${CMAKE_CURRENT_BINARY_DIR}/include/cblas_mangling.h COMPONENT headers
	DESTINATION include)

if(BUILD_TESTS)
  add_subdirectory(testing)
endif(BUILD_TESTS)

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BUILD_EXAMPLES)
