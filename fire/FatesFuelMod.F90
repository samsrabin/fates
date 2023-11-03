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

      procedure :: UpdateLoading

  end type fuel_type

  contains 

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
      
      ! sum across fuel 
      this%total_sum = sum(this%loading(1:nfsc))

      ! update fraction - check for zero total fuel
      if (this%total_sum > 0.0_r8) then 
        this%frac(1:nfsc) = this%loading(1:nfsc)/this%total_sum
      else 
        this%frac(1:nfsc) = 0.0_r8
      end if
        
    end subroutine UpdateLoading

end module FatesFuelMod