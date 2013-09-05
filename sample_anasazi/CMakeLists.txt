# sample program using Anasazi in Trilinos
include("${TRILINOS_PREFIX}/lib/cmake/Anasazi/AnasaziConfig.cmake")

message(STATUS "Display Anasazi_INCLUDE_DIRS=${Anasazi_INCLUDE_DIRS}.")
message(STATUS "Display Anasazi_LIBRARY_DIRS=${Anasazi_LIBRARY_DIRS}.")
message(STATUS "Display Anasazi_LIBRARIES=${Anasazi_LIBRARIES}.")

include_directories(${Anasazi_INCLUDE_DIRS})
link_directories(${Anasazi_LIBRARY_DIRS})

set(TESTS spin_hamiltonian_parallel_LOBPCG spin_hamiltonian_parallel_BlockKrylovSchur_product)

if( DEFINED TRILINOS_PREFIX )
  foreach(name ${TESTS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} ${Anasazi_LIBRARIES})
  endforeach(name)
endif( DEFINED TRILINOS_PREFIX )
