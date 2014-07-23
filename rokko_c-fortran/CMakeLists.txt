set(SOURCES_C wrapper_c.cpp timer_wrapper_c.cpp frank_matrix_wrapper_c.cpp)
add_library(rokko_c ${SOURCES_C})
target_link_libraries(rokko_c rokko)

if(BUILD_Fortran_BINDING)
  set(SOURCES_FORTRAN wrapper_fortran.c timer_wrapper_fortran.c frank_matrix_wrapper_fortran.c rokko_fortran.f90)
endif(BUILD_Fortran_BINDING)
add_library(rokko_fortran ${SOURCES_C} ${SOURCES_FORTRAN})
target_link_libraries(rokko_fortran rokko)

install(TARGETS rokko_c COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(TARGETS rokko_fortran COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(DIRECTORY . COMPONENT headers DESTINATION include/rokko
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
        PATTERN ".git" EXCLUDE PATTERN "CVS" EXCLUDE)

install(DIRECTORY "${PROJECT_BINARY_DIR}/rokko_c-fortran" COMPONENT headers DESTINATION include
        FILES_MATCHING PATTERN "*.mod"
        PATTERN ".git" EXCLUDE PATTERN "CVS" EXCLUDE)
