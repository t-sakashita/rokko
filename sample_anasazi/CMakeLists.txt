# sample program using Anasazi in Trilinos
include("${TRILINOS_PREFIX}/lib/cmake/Anasazi/AnasaziConfig.cmake")

message(STATUS "Display Anasazi_INCLUDE_DIRS=${Anasazi_INCLUDE_DIRS}.")
message(STATUS "Display Anasazi_LIBRARY_DIRS=${Anasazi_LIBRARY_DIRS}.")
message(STATUS "Display Anasazi_LIBRARIES=${Anasazi_LIBRARIES}.")

include_directories(${Anasazi_INCLUDE_DIRS})
link_directories(${Anasazi_LIBRARY_DIRS})

set(TESTS heisenberg_matrix heisenberg_matrix_LOBPCG heisenberg_product_parallel)

if( DEFINED TRILINOS_PREFIX )
  foreach(name ${TESTS})
    set(target_name sample_anasazi_diagonalize_${name})
    message(STATUS "Display target_name=${target_name}.")
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} ${Anasazi_LIBRARIES})
  endforeach(name)
endif( DEFINED TRILINOS_PREFIX )
