module SFFireWeather

  use FatesConstantsMod,   only : r8 => fates_r8

  implicit none
  private 

  type, abstract, public :: fire_weather
    contains 
    procedure(calc_fire_weather_fn), deferred :: Calculate
  end type fire_weather

  abstract interface
    function calc_fire_weather_fn(this)
      import :: fire_weather 
      class(fire_weather), intent(in) :: this 
    end function calc_fire_weather_fn
  end interface 

end module SFFireWeather