program FatesUnitTestSF
  !
  ! DESCRIPTION:
  !		Test the FATES SPITFIRE model
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestSFMod, only : ReadDatmData, WriteFireData
  use FatesUnitTestIOMod, only : logf, OpenFile
  use SFFireWeather,      only : fire_weather
  use SFParamsMod,        only : min_precip_thresh

  implicit none

  ! LOCALS:
  type(fire_weather)    :: fireWeather  ! fire weather object
  real(r8), allocatable :: temp_degC(:) ! daily air temperature [degC]
  real(r8), allocatable :: precip(:)    ! daily precipitation [mm]
  real(r8), allocatable :: rh(:)        ! daily relative humidity [%]
  real(r8), allocatable :: wind(:)      ! daily wind speed [m/s]
  integer,  allocatable :: time(:)      ! time 
  real(r8), allocatable :: NI(:)        ! cumulative nesterov index
  integer               :: n = 365
  integer               :: i

  ! open log file
  logf = OpenFile("log.txt", mode='rw')

  ! allocate arrays
  allocate(temp_degC(n))
  allocate(precip(n))
  allocate(rh(n))
  allocate(wind(n))
  allocate(time(n))
  allocate(NI(n))

  ! read in DATM data
  call ReadDatmData('BONA_datm.nc', temp_degC, precip, rh, wind)

  ! initialize fire weather
  call fireWeather%Init(1)
  
  ! run on time steps
  NI(1) = 0.0
  do i = 1, n
    time(i) = i
    
    ! calculate nesterov index
    if (precip(i) > min_precip_thresh .or. i == 1) then 
      NI(i) = 0.0
    else 
      NI(i) = NI(i-1) + fireWeather%NesterovIndex(temp_degC(i), precip(i), rh(i))
    end if
  end do 

  ! write out data
  call WriteFireData('Fire_unit_out.nc', n, time, temp_degC, precip, rh, NI)

end program FatesUnitTestSF