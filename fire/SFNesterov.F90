module SFNesterov

  use FatesConstantsMod, only : r8 => fates_r8
  use SFFireWeather,     only : fire_weather
  
  implicit none
  private

  type, public, extents(fire_weather) :: NesterovIndex 

    contains 

      procedure, public :: Calculate => calc_nesterov_index

  end type NesterovIndex

  private :: calc_nesterov_index

  contains 

    real(r8) function calc_nesterov_index(this, temp_C, precip, rh)
      !
      !  DESCRIPTION:
      !  Calculates current day's Nesterov Index for a given input values

      use SFParamsMod, only : SF_val_fdi_a, SF_val_fdi_b, min_precip_thresh
      
      !
      ! ARGUMENTS
      class(NesterovIndex), intent(in)   :: this   ! fire weather object
      real(r8),            intent(in)    :: temp_C ! daily averaged temperature [degrees C]
      real(r8),            intent(in)    :: precip ! daily precipitation [mm]
      real(r8),            intent(in)    :: rh     ! daily relative humidity [rh]

      ! LOCALS:
      real(r8) :: yipsolon ! intermediate varable for dewpoint calculation
      real(r8) :: dewpoint ! dewpoint

      if (precip > min_precip_thresh) then ! NI is 0.0 if it rains
        calc_nesterov_index = 0.0_r8
      else 
        ! Calculate dewpoint temperature 
        yipsolon = (SF_val_fdi_a*temp_C)/(SF_val_fdi_b + temp_C) + log(max(1.0_r8, rh)/100.0_r8) 
        dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a - yipsolon) 
        
        ! Nesterov 1968.  Eq 5, Thonicke et al. 2010
        calc_nesterov_index = (temp_C - dewpoint)*temp_C 
        if (calc_nesterov_index < 0.0_r8) calc_nesterov_index = 0.0_r8 ! can't be negative
      endif

    end function calc_nesterov_index

end module SFNesterov
