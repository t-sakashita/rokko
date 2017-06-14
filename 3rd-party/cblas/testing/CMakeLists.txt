add_executable(xscblat1 c_sblat1.f c_sblas1.c)
target_link_libraries(xscblat1 cblas)
add_test(xscblat1, xscblat1)

add_executable(xscblat2 c_sblat2.f c_sblas2.c c_s2chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xscblat2 cblas)

add_executable(xscblat3 c_sblat3.f c_sblas3.c c_s3chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xscblat3 cblas)

add_executable(xdcblat1 c_dblat1.f c_dblas1.c)
target_link_libraries(xdcblat1 cblas)
add_test(xdcblat1, xdcblat1)

add_executable(xdcblat2 c_dblat2.f c_dblas2.c c_d2chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xdcblat2 cblas)

add_executable(xdcblat3 c_dblat3.f c_dblas3.c c_d3chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xdcblat3 cblas)

add_executable(xccblat1 c_cblat1.f c_cblas1.c)
target_link_libraries(xccblat1 cblas)
add_test(xccblat1, xccblat1)

add_executable(xccblat2 c_cblat2.f c_cblas2.c c_c2chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xccblat2 cblas)

add_executable(xccblat3 c_cblat3.f c_cblas3.c c_c3chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xccblat3 cblas)

add_executable(xzcblat1 c_zblat1.f c_zblas1.c)
target_link_libraries(xzcblat1 cblas)
add_test(xzcblat1 xzcblat1)

add_executable(xzcblat2 c_zblat2.f c_zblas2.c c_z2chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xzcblat2 cblas)

add_executable(xzcblat3 c_zblat3.f c_zblas3.c c_z3chke.c auxiliary.c c_xerbla.c)
target_link_libraries(xzcblat3 cblas)
