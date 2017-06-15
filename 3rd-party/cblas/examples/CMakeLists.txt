set(PROG cblas_example1 cblas_example2)
foreach(name ${PROG})
  add_executable(${name} ${name}.c)
  target_link_libraries(${name} cblas)
endforeach(name)

add_test(cblas_example1 cblas_example1)
