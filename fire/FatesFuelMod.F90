module FatesFuelMod
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesLitterMod,    only : ncwd, litter_type

  implicit none 
  private

  ! There are six fuel classes:
  ! 1:4) four CWD_AG pools (twig, s branch, l branch, trunk), 5) dead leaves and 6) live grass
  integer, parameter, public :: nfsc  = ncwd + 2 ! number fuel size classes (4 cwd size classes, leaf litter, and grass)
  integer, parameter, public :: tw_sf = 1        ! array index of twig pool for fire
  integer, parameter, public :: sb_sf = 2        ! array index for small branch pool for fire
  integer, parameter, public :: lb_sf = 3        ! array index of large branch pool for fire
  integer, parameter, public :: tr_sf = 4        ! array index of dead trunk pool for fire
  integer, parameter, public :: dl_sf = 5        ! array index of dead leaf pool for fire (dead grass and dead leaves)
  integer, parameter, public :: lg_sf = 6        ! array index of live grass pool for fire

  type, public :: fuel_type
    
    real(r8) :: loading(nfsc)      ! fuel in individual fuel classes [kgC/m2]
    real(r8) :: total_sum          ! sum across all fuel classes     [kgC/m2]
    real(r8) :: frac(nfsc)         ! fuel fraction across all fuel classes [0-1]
    real(r8) :: moisture(nfsc)     ! fuel moisture across all fuel classes [m3/m3]
    real(r8) :: mean_eff_moisture  ! mean effective fuel moisture [m3/m3]
    real(r8) :: bulk_density       ! average fuel bulk density [kgC/m3]
    real(r8) :: SAV                ! average fuel surface area to volume ratio [/cm]
    real(r8) :: MEF                ! average moisture of extinction [m3/m3]
    real(r8) :: av_moisture        ! average moisture [m3/m3]
    real(r8) :: eff_moisture(nfsc) ! effective moisture, moisture/MEF

    contains 

      procedure :: Init
      procedure :: CalculateLoading
      procedure :: UpdateLoading
      procedure :: FuseFuel
      procedure :: UpdateGeometry
      procedure :: UpdateMoisture

  end type fuel_type

  !=======================================================================================

  contains 

  subroutine Init(this) 
    ! DESCRIPTION:
    !   Initialize fuel class
    
    ! ARGUMENTS:
    class(fuel_type), intent(inout) :: this ! fuel object

    ! just zero everything
    this%loading(:)      = 0.0_r8
    this%total_sum       = 0.0_r8
    this%frac(:)         = 0.0_r8
    this%moisture(:)     = 0.0_r8
    this%bulk_density    = 0.0_r8
    this%SAV             = 0.0_r8
    this%MEF             = 0.0_r8 
    this%av_moisture     = 0.0_r8 
    this%eff_moisture(:) = 0.0_r8 

  end subroutine Init

  !=======================================================================================

  subroutine FuseFuel(this, other_fuel, this_area, other_area)
    ! DESCRIPTION:
    ! Fuse this fuel with other_fuel patch 
    
    ! ARGUMENTS:
    class(fuel_type),  intent(inout) :: this       ! recipient patch's fuel
    class(fuel_type),  intent(in)    :: other_fuel ! donor patch fuel
    real(r8)                         :: this_area  ! area of recipient patch [m2]
    real(r8)                         :: other_area ! area of donor patch [m2]

    ! LOCALS
    integer  :: i            ! looping index
    real(r8) :: inv_sum_area ! inverse sum of patch areas

    inv_sum_area = 1.0_r8/(this_area + other_area)

    do i = 1, nfsc
      this%loading(i) = FuseAttribute(this%loading(i), other_fuel%loading(i),          &
        this_area, other_area)
    end do 

    ! re-sum loading
    call this%UpdateLoading()

    this%bulk_density = FuseAttribute(this%bulk_density, other_fuel%bulk_density,      &
      this_area, other_area)
    this%SAV = FuseAttribute(this%SAV, other_fuel%SAV, this_area, other_area)

  end subroutine FuseFuel

  !=======================================================================================

  real(r8) function FuseAttribute(var, other_var, area, other_area)
    ! DESCRIPTION:
    ! Fuse attributes by area

    ! ARGUMENTS:
    real(r8), intent(in) :: var        ! variable to fuse 
    real(r8), intent(in) :: other_var  ! donor variable
    real(r8), intent(in) :: area       ! area of recipient variable [m2]
    real(r8), intent(in) :: other_area ! area of donor variable [m2]

    ! LOCALS:
    real(r8) :: inv_sum_area ! inverse sum of areas 

    inv_sum_area = 1.0_r8/(area + other_area)

    FuseAttribute = (var*area + other_var*other_area)*inv_sum_area

  end function FuseAttribute

  !=======================================================================================

  subroutine CalculateLoading(this, patch_litter, live_grass)
    ! DESCRIPTION:
    !   Calculate the fuel loading from input litter information
    
    ! ARGUMENTS:
    class(fuel_type),  intent(inout) :: this         ! fuel object
    type(litter_type), intent(in)    :: patch_litter ! patch litter object
    real(r8),          intent(in)    :: live_grass   ! amount of live grass on patch [kg/m2]

    ! coarse woody debris
    this%loading(tw_sf:tr_sf) = patch_litter%ag_cwd(:)

    ! dead leaves 
    this%loading(dl_sf) = sum(patch_litter%leaf_fines(:))

    ! live grass
    this%loading(lg_sf) = live_grass
    
    ! sum up loading
    call this%UpdateLoading()

  end subroutine CalculateLoading

  !=======================================================================================

  subroutine UpdateGeometry(this, bulk_density, sav)
    ! DESCRIPTION:
    !  Update geometry for fuel conditions
    
    ! ARGUMENTS:
    class(fuel_type),  intent(inout) :: this               ! fuel object
    real(r8),          intent(in)    :: bulk_density(nfsc) ! bulk density of fuel [kg/m3]
    real(r8),          intent(in)    :: sav(nfsc)          ! surface area to volume ratio of fuel [/cm]

    ! LOCALS:
    integer :: l ! looping index for litter class

    ! sum up loading - usually we've already done this but do again just in case
    call this%UpdateLoading()

    ! take weighted average across all fuel types (except trunks)
    this%bulk_density = 0.0_r8
    this%SAV = 0.0_r8
    if (this%total_sum > 0.0_r8) then 
      do l = 1, nfsc
        ! don't include trunks
        if (l /= tr_sf) then 
          this%bulk_density = this%bulk_density + this%frac(l)*bulk_density(l)
          this%SAV = this%SAV + this%frac(l)*sav(l)
        end if 
      end do 
    end if 

  end subroutine UpdateGeometry

  !=======================================================================================

  subroutine UpdateMoisture(this, sav)
    ! DESCRIPTION:
    !  Update moisture conditions by weighted average across fuel types
    
    ! ARGUMENTS:
    class(fuel_type),  intent(inout) :: this      ! fuel object
    real(r8),          intent(in)    :: sav(nfsc) ! surface area to volume ratio [/cm]

    ! LOCALS:
    integer  :: l   ! looping index for litter class
    real(r8) :: MEF ! moisture of extinction

    this%MEF = 0.0_r8
    this%av_moisture = 0.0_r8
    if (this%total_sum > 0.0_r8) then
      do l = 1, nfsc
        ! moisture of extinction - moisture at which the fuel cannot burn [m3/m3]
        MEF = MoistureOfExtinction(sav(l))
        
        ! skip trunks for average MEF and moisture
        if (l /= tr_sf) then 
          this%MEF = this%MEF + this%frac(l)*MEF
          this%av_moisture = this%av_moisture + this%frac(l)*this%moisture(l)
        end if 

        this%eff_moisture(l) = this%moisture(l)/MEF
      end do 

    end if 

  end subroutine UpdateMoisture

  !=======================================================================================

  subroutine UpdateLoading(this)
    ! DESCRIPTION:
    !   Sum the loading across all fuel classes and update fuel fraction
    
    ! ARGUMENTS:
    class(fuel_type), intent(inout) :: this ! fuel object

    ! sum across fuel 
    this%total_sum = sum(this%loading(1:nfsc))

    ! update fraction - check for zero total fuel
    if (this%total_sum > 0.0_r8) then 
      this%frac(1:nfsc) = this%loading(1:nfsc)/this%total_sum
    else 
      this%frac(1:nfsc) = 0.0_r8
    end if

  end subroutine UpdateLoading

  !=======================================================================================

  real(r8) function MoistureOfExtinction(sav)
    ! DESCRIPTION:
    !   Calculate moisture of extinction (MEF) [m3/m3]
    !
    
    !   Realistically, MEF depends on compactness of fuel, depth, particle size, wind, slope
    !   Here, we use Eq 27 from Peterson and Ryan (1986) "Modeling Postfire Conifer Mortality for Long-Range Planning"
    !   MEF: pine needles=0.30 (text near Eq 28 Rothermal 1972)
    !      short grass=0.12, tall grass=0.25, chaparral=0.2, closed timber litter=0.30, hardwood litter=0.25
    !   From Thonicke 2010, SAV values in paper give MEF values of:
    !      twigs=0.355; small boles=0.44; large boles=0.525; trunks=0.63; dead/live grass=0.248
    !   Lasslop 2014 Table 1 MEF PFT level:grass=0.2, shrubs=0.3, TropEverGrnTree=0.2, TropDecid Tree=0.3, Extra-trop Tree=0.3

    ! ARGUMENTS:
    real(r8),         intent(in) :: sav ! fuel surface area to volume ratio (/cm)

    if (sav <= 0.0_r8) then 
      MoistureOfExtinction = 0.0_r8
    else
      MoistureOfExtinction = 0.524_r8 - 0.066_r8*log(sav)
    end if

  end function MoistureofExtinction

  !=======================================================================================

end module FatesFuelMod