set(pfunit_sources SFFireWeather.F90)

add_pfunit_ctest(FireWeather
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)

