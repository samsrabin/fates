set(pfunit_sources test_MossFwet.pf test_MossFwetScalerLinear.pf)

add_pfunit_ctest(MossFwet
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  