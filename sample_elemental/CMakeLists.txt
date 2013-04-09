set(TESTS  SymmetricEig_iterator) 
#HermitianEig HermitianEigFromSequential

include_directories(${ELEMENTAL_INCLUDE_DIRS})

foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${ELEMENTAL_LIBRARIES})
  #add_alps_test(${name})
endforeach(name)



