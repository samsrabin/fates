module SFFireWeatherMod

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none
  private 

  type, abstract, public :: fire_weather
    real(r8) :: fire_weather_index   ! fire weather index
    real(r8) :: wind_speed           ! wind speed [m/min]
    real(r8) :: effective_wind_speed ! effective wind speed, corrected for by tree/grass cover [m/min]
    contains 
    procedure(initialize_fire_weather), public, deferred :: Init
    procedure(update_fire_weather),     public, deferred :: Update
    procedure(calculate_fuel_moisture), public, deferred :: CalcFuelMoisture
    procedure,                          public           :: UpdateEffectiveWindSpeed
  end type fire_weather

  abstract interface
    subroutine initialize_fire_weather(this)
      import :: fire_weather 
      class(fire_weather), intent(inout) :: this
    end subroutine initialize_fire_weather
    subroutine update_fire_weather(this, temp_C, precip, rh, wind)
      use FatesConstantsMod, only : r8 => fates_r8
      import :: fire_weather 
      class(fire_weather), intent(inout) :: this
      real(r8),            intent(in)    :: temp_C
      real(r8),            intent(in)    :: precip
      real(r8),            intent(in)    :: rh
      real(r8),            intent(in)    :: wind
    end subroutine update_fire_weather

    real(r8) function calculate_fuel_moisture(this, sav)
      use FatesConstantsMod, only : r8 => fates_r8
      import :: fire_weather
      class(fire_weather), intent(in) :: this
      real(r8),            intent(in) :: sav
    end function calculate_fuel_moisture
  end interface 

  contains

  subroutine UpdateEffectiveWindSpeed(this, tree_fraction, grass_fraction,          &
    bare_fraction)
    !
    !  DESCRIPTION:
    !  Calculates effective wind speed
      
    ! ARGUMENTS
    class(fire_weather), intent(inout) :: this           ! fire weather class
    real(r8),            intent(in)    :: tree_fraction  ! tree fraction [0-1]
    real(r8),            intent(in)    :: grass_fraction ! grass fraction [0-1]
    real(r8),            intent(in)    :: bare_fraction  ! bare ground fraction [0-1]

    this%effective_wind_speed = this%wind_speed*(tree_fraction*0.4_r8 +                  &
      (grass_fraction + bare_fraction)*0.6_r8)

  end subroutine UpdateEffectiveWindSpeed
  
end module SFFireWeatherMod