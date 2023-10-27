module SFFireWeather

  use FatesConstantsMod,   only : r8 => fates_r8

  implicit none
  private 

  type, public :: fire_weather

    integer :: fire_weather_equation ! fire weather equations to use [1=Nesterov; 2=Canadian Fire Weather Index]

    contains 

    procedure :: Init
    procedure :: NesterovIndex

  end type fire_weather

  contains 

    !=====================================================================================

    subroutine Init(this, equation_type)
      !
      !  DESCRIPTION:
      !  Initialize the fire weather class
      !
      ! ARGUMENTS:
      class(fire_weather), intent(inout) :: this          ! fire weather object
      integer,             intent(in)    :: equation_type ! fire weather equations to use [1=Nesterov; 2=Canadian Fire Weather Index]

      this%fire_weather_equation = equation_type

    end subroutine Init

    !=====================================================================================

    subroutine NesterovIndex(this, temp_C, precip, rh, NI)
      !
      !  DESCRIPTION:
      !  Calculates current day's Nesterov Index for a given input values

      use SFParamsMod, only : SF_val_fdi_a, SF_val_fdi_b, min_precip_thresh
      
      !
      ! ARGUMENTS
      class(fire_weather), intent(in)    :: this   ! fire weather object
      real(r8),            intent(in)    :: temp_C ! daily averaged temperature [degrees C]
      real(r8),            intent(in)    :: precip ! daily precipitation [mm]
      real(r8),            intent(in)    :: rh     ! daily relative humidity [rh]
      real(r8),            intent(out)   :: NI   ! daily Nesterov Index [C^2] 

      ! LOCALS:
      real(r8) :: yipsolon ! intermediate varable for dewpoint calculation
      real(r8) :: dewpoint ! dewpoint

      if (precip > min_precip_thresh) then ! NI is 0.0 if it rains
        NI = 0.0_r8
      else 
        
        ! Calculate dewpoint temperature 
        yipsolon = (SF_val_fdi_a*temp_C)/(SF_val_fdi_b + temp_C) + log(max(1.0_r8, rh)/100.0_r8) 
        dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a - yipsolon) 
        
        ! Nesterov 1968.  Eq 5, Thonicke et al. 2010
        NI = (temp_C - dewpoint)*temp_C 
        if (NI < 0.0_r8) NI = 0.0_r8 ! can't be negative
      endif

    end subroutine NesterovIndex

end module SFFireWeather