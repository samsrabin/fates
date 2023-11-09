module FatesFuelMod
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesLitterMod,    only : ncwd, litter_type

  implicit none 
  private

  ! There are six fuel classes:
  ! 1:4) four CWD_AG pools (twig, s branch, l branch, trunk), 5) dead leaves and 6) live grass
  integer, parameter, public :: nfsc  = ncwd + 2 ! number fuel size classes (4 cwd size classes, leaf litter, and grass)
  integer, parameter, public :: tw_sf = 1        ! array index of twig pool for fire
  integer, parameter, public :: lb_sf = 3        ! array index of large branch pool for fire
  integer, parameter, public :: tr_sf = 4        ! array index of dead trunk pool for fire
  integer, parameter, public :: dl_sf = 5        ! array index of dead leaf pool for fire (dead grass and dead leaves)
  integer, parameter, public :: lg_sf = 6        ! array index of live grass pool for fire

  type, public :: fuel_type
    real(r8) :: loading(nfsc) ! fuel in individual fuel classes [kg/m2]
    real(r8) :: total_sum     ! sum across all fuel classes     [kg/m2]
    real(r8) :: frac(nfsc)    ! fuel fraction across all fuel classes [0-1]

    contains 

      procedure :: Init
      procedure :: UpdateLoading
      procedure :: SumLoading
      procedure :: FuseFuel

  end type fuel_type

  !=======================================================================================

  contains 

    subroutine Init(this) 
      ! DESCRIPTION:
      !   Initialize fuel class
      
      ! ARGUMENTS:
      class(fuel_type), intent(inout) :: this ! fuel object

      ! just zero everything
      this%loading(:) = 0.0_r8
      this%total_sum  = 0.0_r8
      this%frac(:)    = 0.0_r8

    end subroutine Init

    !=====================================================================================

    subroutine FuseFuel(this, other_fuel, this_area, other_area)
      ! DESCRIPTION:
      ! Fuse this patch with other_fuel patch 
      
      ! ARGUMENTS:
      class(fuel_type),  intent(inout) :: this       ! recipient patch's fuel
      class(fuel_type),  intent(in)    :: other_fuel ! donor patch fuel
      real(r8)                         :: this_area  ! area of recipient patch [m2]
      real(r8)                         :: other_area ! area of donor patch [m2]

      ! LOCALS
      integer  :: i            ! looping index
      real(r8) :: inv_sum_area ! inverse sum of patch areas

      inv_sum_area = 1.0_r8/(this_area + other_area)

      do i = 1, NFSC
        this%loading(i) = (this%loading(i)*this_area + other_fuel%loading(i)*other_area)*inv_sum_area
      end do 

      ! re-sum loading
      call this%SumLoading()

    end subroutine FuseFuel

    !=====================================================================================

    subroutine UpdateLoading(this, patch_litter, live_grass)
      ! DESCRIPTION:
      !   Update the fuel loading from input litter information
      
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
      call this%SumLoading()

    end subroutine UpdateLoading

    !=====================================================================================

    subroutine SumLoading(this)
      ! DESCRIPTION:
      !   Sum the loading across all fuel classes and update fuel fraction
      
      ! ARGUMENTS:
      class(fuel_type),  intent(inout) :: this         ! fuel object

      ! sum across fuel 
      this%total_sum = sum(this%loading(1:nfsc))

      ! update fraction - check for zero total fuel
      if (this%total_sum > 0.0_r8) then 
        this%frac(1:nfsc) = this%loading(1:nfsc)/this%total_sum
      else 
        this%frac(1:nfsc) = 0.0_r8
      end if

    end subroutine SumLoading

end module FatesFuelMod