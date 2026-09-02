module FatesFuelMod

  use FatesFuelClassesMod, only : fuel_classes
  use FatesInterfaceTypesMod, only : num_fuel_classes
  use FatesInterfaceTypesMod, only : hlm_moss_fuel_moisture_live_intercept
  use FatesInterfaceTypesMod, only : hlm_moss_fuel_moisture_live_slope
  use FatesInterfaceTypesMod, only : hlm_moss_fuel_moisture_dead_intercept
  use FatesInterfaceTypesMod, only : hlm_moss_fuel_moisture_dead_slope
  use FatesConstantsMod,   only : r8 => fates_r8
  use FatesConstantsMod,   only : nearzero
  use SFNesterovMod,       only : nesterov_index
  use SFFireWeatherMod,    only : fire_weather
  use FatesGlobals,        only : fates_log
  use FatesGlobals,        only : endrun => fates_endrun
  use shr_log_mod,         only : errMsg => shr_log_errMsg

  implicit none
  private

  type, public :: fuel_type
    
    real(r8), allocatable :: loading(:)            ! fuel loading of each fuel class [kgC/m2]
    real(r8), allocatable :: effective_moisture(:) ! fuel effective moisture all fuel class (moisture/MEF) [m3/m3]
    real(r8), allocatable :: frac_loading(:)       ! fractional loading of all fuel classes [0-1]
    real(r8), allocatable :: frac_burnt(:)         ! fraction of litter burnt by fire [0-1]
    real(r8) :: non_trunk_loading                    ! total fuel loading excluding trunks [kgC/m2]
    real(r8) :: average_moisture_notrunks            ! weighted average of fuel moisture across non-trunk fuel classes [m3/m3]
    real(r8) :: bulk_density_notrunks                ! weighted average of bulk density across non-trunk fuel classes [kg/m3]
    real(r8) :: SAV_notrunks                         ! weighted average of surface area to volume ratio across non-trunk fuel classes [/cm]
    real(r8) :: MEF_notrunks                         ! weighted average of moisture of extinction across non-trunk fuel classes [m3/m3]

    contains
      
      procedure :: Init
      procedure :: Fuse
      procedure :: UpdateLoading
      procedure :: SumLoading
      procedure :: CalculateFractionalLoading
      procedure :: UpdateFuelMoisture
      procedure :: AverageBulkDensity_NoTrunks
      procedure :: AverageSAV_NoTrunks
      procedure :: CalculateFuelBurnt
      procedure :: CalculateResidenceTime
      procedure :: Deallocate

  end type fuel_type

  ! CONSTANTS:
  ! Public so that unit tests can assert against it rather than restating its value.
  real(r8), parameter, public :: max_grass_frac = 0.8_r8 ! maximum fraction burnt for live grass fuels

  ! Public for the same reason: unit tests need the moisture of extinction that
  ! UpdateFuelMoisture normalizes by, and reimplementing it in the test would let the two
  ! drift apart silently.
  public :: MoistureOfExtinction

  contains 

    subroutine Init(this)
      ! DESCRIPTION:
      !   Initialize fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      ! Allocate
      allocate(this%loading(1:num_fuel_classes))
      allocate(this%frac_loading(num_fuel_classes))
      allocate(this%frac_burnt(num_fuel_classes))
      allocate(this%effective_moisture(num_fuel_classes))
      
      ! just zero everything
      this%loading(1:num_fuel_classes) = 0.0_r8
      this%frac_loading(1:num_fuel_classes) = 0.0_r8
      this%frac_burnt(1:num_fuel_classes) = 0.0_r8  
      this%effective_moisture(1:num_fuel_classes) = 0.0_r8
      this%non_trunk_loading = 0.0_r8
      this%average_moisture_notrunks = 0.0_r8 
      this%bulk_density_notrunks = 0.0_r8
      this%SAV_notrunks = 0.0_r8
      this%MEF_notrunks = 0.0_r8 

    end subroutine Init 
    
    !-------------------------------------------------------------------------------------
    
    subroutine Fuse(this, self_area, donor_area, donor_fuel)
      ! DESCRIPTION:
      !   Fuse attributes of this object with another

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this       ! fuel class
      real(r8),         intent(in)    :: self_area  ! area of this fuel class's patch [m2]
      real(r8),         intent(in)    :: donor_area ! area of donor fuel class's patch [m2]
      type(fuel_type),  intent(in)    :: donor_fuel ! donor fuel class
      
      ! LOCALS:
      integer  :: i            ! looping index
      real(r8) :: self_weight  ! weighting of the receiving fuel class
      real(r8) :: donor_weight ! weighting of the donor fuel class
      
      self_weight = self_area/(donor_area + self_area)
      donor_weight = 1.0_r8 - self_weight
      
      do i = 1, num_fuel_classes 
        this%loading(i) = this%loading(i)*self_weight +                                  &
          donor_fuel%loading(i)*donor_weight
        this%frac_loading(i) = this%frac_loading(i)*self_weight +                        &
          donor_fuel%frac_loading(i)*donor_weight
        this%frac_burnt(i) = this%frac_burnt(i)*self_weight +                            &
          donor_fuel%frac_burnt(i)*donor_weight
        this%effective_moisture(i) = this%effective_moisture(i)*self_weight +            &
          donor_fuel%effective_moisture(i)*donor_weight
      end do 
      
      this%non_trunk_loading = this%non_trunk_loading*self_weight +                      &
        donor_fuel%non_trunk_loading*donor_weight
      this%average_moisture_notrunks = this%average_moisture_notrunks*self_weight +                        &
        donor_fuel%average_moisture_notrunks*donor_weight
      this%bulk_density_notrunks = this%bulk_density_notrunks*self_weight +              &
        donor_fuel%bulk_density_notrunks*donor_weight      
      this%SAV_notrunks = this%SAV_notrunks*self_weight + donor_fuel%SAV_notrunks*donor_weight
      this%MEF_notrunks = this%MEF_notrunks*self_weight + donor_fuel%MEF_notrunks*donor_weight    
      
    end subroutine Fuse

    !-------------------------------------------------------------------------------------

    subroutine UpdateLoading(this, leaf_litter, twig_litter, small_branch_litter,     &
        large_branch_litter, trunk_litter, live_grass, live_moss, dead_moss)
      ! DESCRIPTION:
      !   Updates loading for each fuel type

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                ! fuel class
      real(r8),         intent(in)    :: leaf_litter         ! input leaf litter [kgC/m2]
      real(r8),         intent(in)    :: twig_litter         ! input twig litter [kgC/m2]
      real(r8),         intent(in)    :: small_branch_litter ! input small branch litter [kgC/m2]
      real(r8),         intent(in)    :: large_branch_litter ! input leaf litter [kgC/m2]
      real(r8),         intent(in)    :: trunk_litter        ! input leaf litter [kgC/m2]
      real(r8),         intent(in)    :: live_grass          ! input live grass [kgC/m2]
      real(r8),         intent(in)    :: live_moss           ! input live moss [kgC/m2]
      real(r8),         intent(in)    :: dead_moss           ! input dead moss [kgC/m2]

      this%loading(fuel_classes%dead_leaves()) = leaf_litter
      this%loading(fuel_classes%twigs()) = twig_litter
      this%loading(fuel_classes%small_branches()) = small_branch_litter
      this%loading(fuel_classes%large_branches()) = large_branch_litter
      this%loading(fuel_classes%live_grass()) = live_grass
      this%loading(fuel_classes%trunks()) = trunk_litter

      ! Special handling for moss because its fuel classes are only allocated if hlm_use_moss true.
      ! Dead-moss loading can be nonzero even when live-moss loading is zero, so both are set
      ! together whenever the moss fuel classes are present.
      if (fuel_classes%moss_classes_present()) then
         this%loading(fuel_classes%live_moss()) = live_moss
         this%loading(fuel_classes%dead_moss()) = dead_moss
      else if (live_moss > 0._r8 .or. dead_moss > 0._r8) then
         write(fates_log(), *) 'Moss biomass > 0. live_moss: ', live_moss, ' dead_moss: ', dead_moss
         write(fates_log(), *) 'but moss fuel classes not present.'
         call endrun(msg=errMsg( __FILE__, __LINE__))
      end if

    end subroutine UpdateLoading

    !-------------------------------------------------------------------------------------

    subroutine SumLoading(this)
      ! DESCRIPTION:
      !   Sums up the loading - excludes trunks
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread 
      !    Rothermel, 1972 (USDA FS GTR INT-115) 
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! LOCALS:
      integer :: i ! looping index
      
      this%non_trunk_loading = 0.0_r8
      do i = 1, num_fuel_classes
        if (i /= fuel_classes%trunks()) then 
          this%non_trunk_loading = this%non_trunk_loading + this%loading(i)
        end if
      end do

    end subroutine SumLoading

    !-------------------------------------------------------------------------------------
    
    subroutine CalculateFractionalLoading(this)
      ! DESCRIPTION:
      !   Calculates fractional loading for fuel

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class
      
      ! LOCALS:
      integer :: i ! looping index

      ! sum up loading just in case
      call this%SumLoading()
      
      if (this%non_trunk_loading > nearzero) then
        do i = 1, num_fuel_classes 
          if (i /= fuel_classes%trunks()) then 
            this%frac_loading(i) = this%loading(i)/this%non_trunk_loading
          else 
            this%frac_loading(i) = 0.0_r8
          end if 
        end do 
      else 
        this%frac_loading(1:num_fuel_classes) = 0.0_r8
        this%non_trunk_loading = 0.0_r8
      end if 

    end subroutine CalculateFractionalLoading

    !-------------------------------------------------------------------------------------
    
    subroutine UpdateFuelMoisture(this, sav_fuel, drying_ratio, fwet_moss, fireWeatherClass)
      ! DESCRIPTION:
      !   Updates fuel moisture depending on what fire weather class is in use
      
      ! ARGUMENTS:
      class(fuel_type),    intent(inout) :: this                       ! fuel class
      real(r8),            intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      real(r8),            intent(in)    :: drying_ratio               ! drying ratio
      real(r8),            intent(in)    :: fwet_moss                  ! moss wetness proxy; ignored if moss fuel classes absent [0-1]
      class(fire_weather), intent(in)    :: fireWeatherClass           ! fireWeatherClass
      
      real(r8) :: moisture(num_fuel_classes)               ! fuel moisture [m3/m3]
      real(r8) :: moisture_of_extinction(num_fuel_classes) ! fuel moisture of extinction [m3/m3]
      integer  :: i                                        ! looping index
 
      if (this%non_trunk_loading + this%loading(fuel_classes%trunks()) > nearzero) then 
        ! calculate fuel moisture [m3/m3] for each fuel class depending on what
        ! fire weather class is in use
        select type (fireWeatherClass)
          class is (nesterov_index)
            call CalculateFuelMoistureNesterov(sav_fuel, drying_ratio,                   &
              fireWeatherClass%fire_weather_index, moisture)
          class default 
            write(fates_log(), *) 'Unknown fire weather class selected.'
            write(fates_log(), *) 'Choose a different fire weather class or upate this subroutine.'
            call endrun(msg=errMsg( __FILE__, __LINE__))
        end select
        
        ! Moss fuel moisture is diagnosed from the moss wetness proxy rather than from the
        ! fire weather index, so the two moss classes are overwritten here. This sits ahead
        ! of the loop below so that moss gets moisture-of-extinction-normalized like every
        ! other fuel class.
        !
        ! The intercepts are constrained, and the constraint is easy to violate by
        ! accident. Every other fuel class takes its moisture from
        ! CalculateFuelMoistureNesterov, which decays toward zero without a floor and so
        ! crosses its moisture of extinction eventually, however wet it starts. This linear
        ! map is instead floored at its intercept. An intercept at or above the class's MEF
        ! therefore makes that class permanently non-flammable rather than merely hard to
        ! ignite: effective_moisture never falls below 1, and CalculateFuelBurnt's "very
        ! wet litter" branch zeroes frac_burnt at every value of fwet_moss.
        !
        ! MEF is MEF_a - MEF_b*log(sav_fuel) (MoistureOfExtinction, below), so for moss it
        ! is fixed by the moss entries of fates_fire_SAV on the FATES parameter file and by
        ! nothing else, rising as SAV falls: at the shipped SAV of 66 /cm it is 0.2475
        ! m3/m3. A moss class burns while fwet_moss is under
        ! (MEF - hlm_moss_fuel_moisture_*_intercept)/hlm_moss_fuel_moisture_*_slope, both
        ! set by the host; at the shipped intercept of 0 and slope of 0.7 that puts both
        ! moss classes' extinction at fwet_moss = 0.354.
        if (fuel_classes%moss_classes_present()) then
          moisture(fuel_classes%live_moss()) = max(0.0_r8,                              &
            hlm_moss_fuel_moisture_live_intercept + hlm_moss_fuel_moisture_live_slope*fwet_moss)
          moisture(fuel_classes%dead_moss()) = max(0.0_r8,                              &
            hlm_moss_fuel_moisture_dead_intercept + hlm_moss_fuel_moisture_dead_slope*fwet_moss)
        end if

        this%average_moisture_notrunks = 0.0_r8
        this%MEF_notrunks = 0.0_r8
        do i = 1, num_fuel_classes
          ! calculate moisture of extinction and fuel effective moisture
          moisture_of_extinction(i) = MoistureOfExtinction(sav_fuel(i))
          this%effective_moisture(i) = moisture(i)/moisture_of_extinction(i)
          
          ! average fuel moisture  and MEF across all fuel types except trunks [m3/m3]
          if (i /= fuel_classes%trunks()) then 
            this%average_moisture_notrunks = this%average_moisture_notrunks + this%frac_loading(i)*moisture(i)
            this%MEF_notrunks = this%MEF_notrunks + this%frac_loading(i)*moisture_of_extinction(i)
          end if 
        end do

      else 
        this%effective_moisture(1:num_fuel_classes) = 0.0_r8
        this%average_moisture_notrunks = 0.0_r8
        this%MEF_notrunks = 0.0_r8
      end if
      
    end subroutine UpdateFuelMoisture
    
    !-------------------------------------------------------------------------------------
        
    subroutine CalculateFuelMoistureNesterov(sav_fuel, drying_ratio, NI, moisture)
      !
      ! DESCRIPTION:
      !   Updates fuel moisture

      ! ARGUMENTS:
      real(r8), intent(in)  :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      real(r8), intent(in)  :: drying_ratio               ! drying ratio
      real(r8), intent(in)  :: NI                         ! Nesterov Index
      real(r8), intent(out) :: moisture(num_fuel_classes) ! moisture of litter [m3/m3]
      
      ! LOCALS
      integer  :: i         ! looping index
      real(r8) :: alpha_FMC ! intermediate variable for calculating fuel moisture
      
      do i = 1, num_fuel_classes
        if (i == fuel_classes%live_grass()) then 
          ! live grass moisture is a function of SAV and changes via Nesterov Index
          ! along the same relationship as the 1 hour fuels
          ! live grass has same SAV as dead grass, but retains more moisture with this calculation
          alpha_FMC = sav_fuel(fuel_classes%twigs())/drying_ratio
        else
          alpha_FMC = sav_fuel(i)/drying_ratio
        end if
        moisture(i) = exp(-1.0_r8*alpha_FMC*NI)
      end do
      
    end subroutine CalculateFuelMoistureNesterov
    
    !-------------------------------------------------------------------------------------
    
    real(r8) function MoistureOfExtinction(sav)
      !
      ! DESCRIPTION:
      !   Calculates moisture of extinction based on input surface area to volume ratio
    
      ! MEF (moisure of extinction) depends on compactness of fuel, depth, particle size, 
      !  wind, and slope
      ! Equation here is Eq. 27 from Peterson and Ryan (1986) "Modeling Postfire Conifer 
      ! Mortality for Long-Range Planning"
      !
      ! Example MEFs:
      ! pine needles = 0.30 (Rothermel 1972)
      ! short grass = 0.12 (Rothermel 1983; Gen. Tech. Rep. INT-143; Table II-1)
      ! tall grass = 0.24 (Rothermel 1983)
      ! chaparral = 0.20 (Rothermel 1983)
      ! closed timber litter = 0.30 (Rothermel 1983)
      ! hardwood litter = 0.25 (Rothermel 1983)
      ! grass = 0.2 (Lasslop 2014; Table 1)
      ! shrubs = 0.3 (Lasslop 2014; Table 1)
      ! tropical evergreen trees = 0.2 (Lasslop 2014; Table 1)
      ! tropical deciduous trees = 0.3 (Lasslop 2014; Table 1)
      ! extratropical trees = 0.3 (Lasslop 2014; Table 1)
      !
      ! SAV values from Thonicke 2010 give: 
      ! twigs = 0.355, small branches = 0.44, large branches = 0.525, trunks = 0.63
      ! dead leaves = 0.248, live grass = 0.248
      !
    
      ! ARGUMENTS:
      real(r8), intent(in) :: sav ! fuel surface area to volume ratio [/cm]
      
      ! CONSTANTS:
      real(r8), parameter :: MEF_a = 0.524_r8
      real(r8), parameter :: MEF_b = 0.066_r8
      
      if (sav <= 0.0_r8) then
        write(fates_log(), *) 'SAV cannot be negative - SAV'
        call endrun(msg=errMsg(__FILE__, __LINE__))
      else
        MoistureOfExtinction = MEF_a - MEF_b*log(sav)
      end if
    
    end function MoistureOfExtinction
    
    !-------------------------------------------------------------------------------------
    
    subroutine AverageBulkDensity_NoTrunks(this, bulk_density)
      ! DESCRIPTION:
      !   Calculates average bulk density (not including trunks)
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread 
      !    Rothermel, 1972 (USDA FS GTR INT-115) 
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type),   intent(inout) :: this                           ! fuel class
      real(r8),           intent(in)    :: bulk_density(num_fuel_classes) ! bulk density of all fuel types [kg/m3]
      
      ! LOCALS:
      integer :: i ! looping index
      
      if (this%non_trunk_loading > nearzero) then
        this%bulk_density_notrunks = 0.0_r8
        do i = 1, num_fuel_classes               
          ! average bulk density across all fuel types except trunks 
          if (i /= fuel_classes%trunks()) then 
            this%bulk_density_notrunks = this%bulk_density_notrunks + this%frac_loading(i)*bulk_density(i)
          end if 
        end do
      else 
        this%bulk_density_notrunks = sum(bulk_density(1:num_fuel_classes))/num_fuel_classes
      end if
      
    end subroutine AverageBulkDensity_NoTrunks
    
    !-------------------------------------------------------------------------------------
    
    subroutine AverageSAV_NoTrunks(this, sav_fuel)
      ! DESCRIPTION:
      !   Calculates average surface area to volume ratio (not including trunks)
      !
      !   Only the 1-h, 10-h and 100-h fuel classes influence fire spread 
      !    Rothermel, 1972 (USDA FS GTR INT-115) 
      !    Wilson, 1982 (UTINT-289)
      !    Pyne et al., 1996 (Introduction to wildland fire)

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                       ! fuel class
      real(r8),         intent(in)    :: sav_fuel(num_fuel_classes) ! surface area to volume ratio of all fuel types [/cm]
      
      ! LOCALS:
      integer :: i ! looping index
      
      if (this%non_trunk_loading > nearzero) then
        this%SAV_notrunks = 0.0_r8
        do i = 1, num_fuel_classes               
          ! average bulk density across all fuel types except trunks 
          if (i /= fuel_classes%trunks()) then 
            this%SAV_notrunks = this%SAV_notrunks + this%frac_loading(i)*sav_fuel(i)
          end if 
        end do
      else 
        this%SAV_notrunks = sum(sav_fuel(1:num_fuel_classes))/num_fuel_classes 
      end if
    
    end subroutine AverageSAV_NoTrunks
    
  !---------------------------------------------------------------------------------------
    
    subroutine CalculateFuelBurnt(this, fuel_consumed)
      ! DESCRIPTION:
      !   Calculates the fraction and total amount of fuel burnt
      !
      
      use SFParamsMod, only : SF_val_mid_moisture, SF_val_mid_moisture_Coeff
      use SFParamsMod, only : SF_val_mid_moisture_Slope, SF_val_min_moisture
      use SFParamsMod, only : SF_val_low_moisture_Coeff, SF_val_low_moisture_Slope
      use SFParamsMod, only : SF_val_miner_total
      use FatesInterfaceTypesMod, only : hlm_moss_max_burn_frac

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this                            ! fuel class
      real(r8),         intent(out)   :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]
      
      ! LOCALS:
      real(r8) :: rel_moisture  ! relative moisture of fuel (moist/moisture of extinction) [unitless]
      integer  :: i             ! looping index
      integer  :: live_moss_i   ! index of the live moss fuel class; 0 if the moss classes are absent
      
      ! Resolve the live moss class index once. It is looked up here rather than inside the
      ! loop because fuel_classes%live_moss() aborts when the moss fuel classes are absent,
      ! and Fortran does not guarantee short-circuit evaluation of .and.
      live_moss_i = 0
      if (fuel_classes%moss_classes_present()) live_moss_i = fuel_classes%live_moss()

      this%frac_burnt(:) = 1.0_r8
        
      ! Calculate fraction of litter is burnt for all classes. 
      ! Equation B1 in Thonicke et al. 2010
      do i = 1, num_fuel_classes        
        
        rel_moisture = this%effective_moisture(i)                  
        
        if (rel_moisture <= SF_val_min_moisture(i)) then
          ! very dry litter
          this%frac_burnt(i) = 1.0_r8 
        else if (rel_moisture > SF_val_min_moisture(i) .and. rel_moisture <= SF_val_mid_moisture(i)) then
          ! low to medium moisture
          this%frac_burnt(i) = max(0.0_r8, min(1.0_r8, SF_val_low_moisture_Coeff(i) - &
            SF_val_low_moisture_Slope(i)*rel_moisture))
        else if (rel_moisture > SF_val_mid_moisture(i) .and. rel_moisture <= 1.0_r8) then
          ! medium to high moisture
          this%frac_burnt(i) = max(0.0_r8, min(1.0_r8, SF_val_mid_moisture_Coeff(i) - &
            SF_val_mid_moisture_Slope(i)*rel_moisture))
        else 
          ! very wet litter
          this%frac_burnt(i) = 0.0_r8  
        endif
        
        ! we can't ever kill all of the grass
        if (i == fuel_classes%live_grass()) then
          this%frac_burnt(i) = min(max_grass_frac, this%frac_burnt(i))
        else if (i == live_moss_i) then
          ! Moss does not inherit the grass cap above, which exists because we can't ever
          ! kill all of the grass; moss carries no such assumption, as a moss mat can burn
          ! off completely. Live moss gets its own limit instead, which defaults to 1.0, at
          ! which this min is a no-op; the knob exists so that the cap can be tightened
          ! during tuning without a code change. Note the capped value is still scaled by
          ! (1 - SF_val_miner_total) below, exactly as grass's cap is, so a cap dialed in
          ! during tuning is not the fraction that finally reaches leaf_burn_frac.
          this%frac_burnt(i) = min(hlm_moss_max_burn_frac, this%frac_burnt(i))
        end if
        
        ! reduce fraction burnt based on mineral content
        this%frac_burnt(i) = this%frac_burnt(i)*(1.0_r8 - SF_val_miner_total)
        
        ! calculate fuel consumed
        fuel_consumed(i) = this%frac_burnt(i)*this%loading(i)
      end do

    end subroutine CalculateFuelBurnt
    
    !-------------------------------------------------------------------------------------
    
    subroutine CalculateResidenceTime(this, tau_l)
      !
      !  DESCRIPTION:
      !  Calculates fire residence time, duration of lethal bole heating [min]
      !  This is used for determining cambial kill of woody cohorts
      !
      !  From Peterson & Ryan (1986)
      !
      
      ! ARGUMENTS:
      class(fuel_type), intent(in)  :: this  ! fuel class
      real(r8),         intent(out) :: tau_l ! duration of lethal bole heating [min]
      
      ! LOCALS:
      integer :: i ! looping index
      
      tau_l = 0.0_r8 
      do i = 1, num_fuel_classes
        if (i /= fuel_classes%trunks()) then 
          ! don't include 1000-hr fuels
          ! convert loading from kgC/m2 to g/cm2
          tau_l = tau_l + 39.4_r8*(this%frac_loading(i)*this%non_trunk_loading/0.45_r8/10.0_r8)* &
            (1.0_r8 - ((1.0_r8 - this%frac_burnt(i))**0.5_r8))
        end if 
      end do
      
      ! cap the residence time to 8mins, as suggested by literature survey by P&R (1986)
      tau_l = min(8.0_r8, tau_l) 

    end subroutine CalculateResidenceTime

    subroutine Deallocate(this)
      ! DESCRIPTION:
      !   Deallocate fuel class

      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel class

      if (allocated(this%loading)) then
         deallocate(this%loading)
      end if
      if (allocated(this%frac_loading)) then
         deallocate(this%frac_loading)
      end if
      if (allocated(this%frac_burnt)) then
         deallocate(this%frac_burnt)
      end if
      if (allocated(this%effective_moisture)) then
         deallocate(this%effective_moisture)
      end if

    end subroutine Deallocate
    
end module FatesFuelMod
