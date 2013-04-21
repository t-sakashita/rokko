cmake_minimum_required(VERSION 2.8)

set(TESTS  diagonalize_elemental SymmetricEig)
# HermitianEigFromSequential Constructors)
#SymmetricEig_iterator)
#HermitianEig HermitianEigFromSequential

include_directories(${ELEMENTAL_INCLUDE_DIRS})

foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${ELEMENTAL_LIBRARIES}  "-mkl=parallel -lifcore")
  #add_alps_test(${name})
endforeach(name)



