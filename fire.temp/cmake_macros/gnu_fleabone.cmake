set(NETCDF_C_PATH "/usr/local/Cellar/netcdf/4.9.2_1")
set(NETCDF_FORTRAN_PATH "/usr/local/Cellar/netcdf-fortran/4.6.1")

# These LDFLAGS provide lapack and blas support on a Mac. This may require installation of
# the Apple Developer Tools.
string(APPEND LDFLAGS " -framework Accelerate")

# This is only needed for the unit test build of the share and/or cpl7 code
string(APPEND LDFLAGS " -L/usr/local/Cellar/netcdf-fortran/4.6.1/lib -lnetcdff")
string(APPEND LDFLAGS " -L/usr/local/Cellar/netcdf/4.9.2_1/lib -lnetcdf")

# Trying to produce a backtrace leads to a hang, so don't even try
string(APPEND CFLAGS " -fno-backtrace -Wno-error=implicit-function-declaration")
# (2023-04-19) The addition of -fallow-argument-mismatch was needed to get some unit tests
# to build. I'm not sure why that's needed given that it appears to already be added by
# the standard gnu.cmake.
string(APPEND FFLAGS " -fno-backtrace -fallow-argument-mismatch")

# This is needed to run the Fortran unit tests;
# this isn't needed to build and run CESM.
if (MPILIB STREQUAL mpi-serial AND NOT compile_threaded)
  set(PFUNIT_PATH "/Users/afoster/pFUnit/build/installed")
endif()

# Most of the following paths aren't necessary on my machine because I have my PATH set so
# that the right compilers are picked up by default. But it doesn't hurt to be explicit.
set(SFC "/usr/local/bin/gfortran")
set(SCC "/usr/bin/gcc")
set(SCXX "/usr/bin/g++")
set(MPIFC "/usr/local/bin/mpif90")
set(MPICC "/usr/local/bin/mpicc")
set(MPICXX "/usr/local/bin/mpic++")
