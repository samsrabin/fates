# This file helps to build unit test programs to test FATES, but is not used in 
# production runs

cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)

project(FATES_tests Fortran C)

include(CIME_utils)

set(HLM_ROOT "../../")

# Add source directories from other share code (csm_share, etc.)
add_subdirectory(${HLM_ROOT}/share/src csm_share)
add_subdirectory(${HLM_ROOT}/share/unit_test_stubs/util csm_share_stubs)

# Add FATES source directories
add_subdirectory(${HLM_ROOT}/src/fates/main fates_main)
add_subdirectory(${HLM_ROOT}/src/fates/biogeochem fates_biogeochem)
add_subdirectory(${HLM_ROOT}/src/fates/fire fates_fire)
add_subdirectory(${HLM_ROOT}/src/fates/radiation fates_radiation)

# Remove shr_mpi_mod from share_sources.
# This is needed because we want to use the mock shr_mpi_mod in place of the real one
#
# TODO: this should be moved into a general-purpose function in Sourcelist_utils.
# Then this block of code could be replaced with a single call, like:
# remove_source_file(${share_sources} "shr_mpi_mod.F90")
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_mpi_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# Remove shr_cal_mod from share_sources.
#
# shr_cal_mod depends on ESMF (or the lightweight esmf wrf timemgr, at
# least). Since CTSM doesn't currently use shr_cal_mod, we're avoiding
# the extra overhead of including esmf_wrf_timemgr sources in this
# build.
#
# TODO: like above, this should be moved into a general-purpose function
# in Sourcelist_utils.  Then this block of code could be replaced with a
# single call, like: remove_source_file(${share_sources}
# "shr_cal_mod.F90")
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_cal_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# Build libraries containing stuff needed for the unit tests.
# Eventually, these add_library calls should probably be distributed into the correct location, rather than being in this top-level CMakeLists.txt file.
add_library(csm_share ${share_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
add_library(fates ${fates_sources})
add_dependencies(fates csm_share)

# We need to look for header files here, in order to pick up shr_assert.h
include_directories(${HLM_ROOT}/share/include)

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the test directories
# Note: it's possible that these could be added by each source directory that
# has tests in it. However, it appears that the order needs to be done
# carefully: for example, include_directories and link_directories needs to be
# done before adding the tests themselves.
add_subdirectory(${HLM_ROOT}/src/fates/fire/test fates_fire_test)