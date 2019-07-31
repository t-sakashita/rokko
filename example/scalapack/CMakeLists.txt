if(BUILD_Fortran_BINDING)
  set(PROG pdsyev_f pdsyevd_f)
  foreach(name ${PROG})
    set(target_name example_lapack_${name})
    add_executable(${target_name} ${name}.f90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
    add_test(${name} ${name})
  endforeach(name)
  unset(PROG)
endif(BUILD_Fortran_BINDING)