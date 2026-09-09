set(pfunit_sources test_RootProfile.pf)

add_pfunit_ctest(RootProfile
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  