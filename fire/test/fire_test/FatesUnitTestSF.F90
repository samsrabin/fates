program FatesUnitTestSF
  !
  ! DESCRIPTION:
  !		Test the FATES SPITFIRE model
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesConstantsMod,  only : fates_unset_int
  use FatesUnitTestSFMod, only : ReadDatmData, WriteFireData, SetUpSite, SetUpParams
  use FatesUnitTestIOMod, only : logf, OpenFile
  use EDTypesMod,         only : ed_site_type
  use FatesPatchMod,      only : fates_patch_type

  use FatesLitterMod,     only : litter_type
  use SFMainMod,          only : UpdateFuelCharacteristics
  use FatesFuelMod,       only : nfsc
  use SFParamsMod,        only : SF_val_fdi_a

  implicit none

  ! LOCALS:
  type(ed_site_type),     pointer     :: site         ! site object
  real(r8),               allocatable :: temp_degC(:) ! daily air temperature [degC]
  real(r8),               allocatable :: precip(:)    ! daily precipitation [mm]
  real(r8),               allocatable :: rh(:)        ! daily relative humidity [%]
  real(r8),               allocatable :: wind(:)      ! daily wind speed [m/s]
  integer,                allocatable :: time(:)      ! time 
  real(r8),               allocatable :: NI(:)        ! cumulative nesterov index
  real(r8),               allocatable :: loading(:)   ! fuel loading [kg/m2]
  real(r8),               allocatable :: moisture(:,:) ! fuel moisture [m3/m3]
  integer                             :: n = 365      ! number of days to run
  integer                             :: i            ! looping indices

  ! open log file
  logf = OpenFile("log.txt", mode='rw')

  ! set up the site and patch
  allocate(site)
  call SetUpSite(site)

  call SetUpParams()

  ! allocate arrays
  allocate(temp_degC(n))
  allocate(precip(n))
  allocate(rh(n))
  allocate(wind(n))
  allocate(time(n))
  allocate(NI(n))
  allocate(loading(n))
  allocate(moisture(n, nfsc))

  ! read in DATM data
  call ReadDatmData('BONA_datm.nc', temp_degC, precip, rh, wind)

  ! run on time steps
  do i = 1, n
    
    time(i) = i
    
    call site%fireWeather%Update(temp_degC(i), precip(i), rh(i))
    NI(i) = site%fireWeather%fire_weather_index
    
    call UpdateFuelCharacteristics(site)
    loading(i) = site%youngest_patch%fuel%total_sum
    moisture(i, 1:nfsc) = site%youngest_patch%fuel%moisture(1:nfsc)
    
  end do 

  ! write out data
  call WriteFireData('Fire_unit_out.nc', n, time, temp_degC, precip, rh, NI, loading,    &
    moisture)

end program FatesUnitTestSF