if(BUILD_EIGEN_SX)
  set(SOURCES communication_sx.f pdlaed2x.f pdlaed3x.f lapack_eigen.f eigen_pdlaedz.f eigen_pdlaed0.f eigen_pdsxedc.f eigen_pdlaed1.f eigen_pbk.f eigen_pbk_main1.f eigen_pbk_main2.f eigen_pbk_step.f machine_epsilon.f matrix_adjust_sx.f eigen_prd.f eigen_prd_2update.f eigen_prd_Au.f eigen_prd_v.f eigen_prd_u.f eigen_prd_2update_local.f eigen_prd_2update_v.f eigen_prd_d.f eigen_prd_init_final.f cstab.f get_delta.c eigen_dcx.f eigen_sx.f eigen_sx_wrap.f)
  add_library(eigen_sx ${SOURCES})
endif(BUILD_EIGEN_SX)
