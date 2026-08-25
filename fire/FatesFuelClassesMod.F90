module FatesFuelClassesMod

  use FatesLitterMod, only : ncwd
  use FatesInterfaceTypesMod, only : num_fuel_classes
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun
  use shr_log_mod, only : errMsg => shr_log_errMsg

  implicit none
  private

  type :: fuel_classes_type
    ! There are eight fuel classes:
    ! 1) twigs, 2) small branches, 3) large branches 4) trunks
    ! 5) dead leaves, 6) live grass, 7) live moss, 8) dead moss
    integer, private :: twigs_i = 1          ! array index for twigs pool
    integer, private :: small_branches_i = 2 ! array index for small branches pool
    integer, private :: large_branches_i = 3 ! array index for large branches pool
    integer, private :: dead_leaves_i = 5    ! array index for dead leaves pool
    integer, private :: live_grass_i = 6     ! array index for live grass pool
    integer, private :: trunks_i = 4         ! array index for trunks pool
    integer, private :: live_moss_i = 7      ! array index for live moss pool
    integer, private :: dead_moss_i = 8      ! array index for dead moss pool

    contains 

    procedure :: twigs, small_branches, large_branches, trunks
    procedure :: dead_leaves, live_grass
    procedure :: live_moss, dead_moss
    procedure :: moss_classes_present

  end type fuel_classes_type

  ! actual type we can pass around 
  type(fuel_classes_type), public :: fuel_classes

  contains

  integer function twigs(this)
    class(fuel_classes_type), intent(in) :: this 
    twigs = this%twigs_i
  end function twigs 

  integer function small_branches(this)
    class(fuel_classes_type), intent(in) :: this 
    small_branches = this%small_branches_i
  end function small_branches 

  integer function large_branches(this)
    class(fuel_classes_type), intent(in) :: this 
    large_branches = this%large_branches_i
  end function large_branches 

  integer function trunks(this)
    class(fuel_classes_type), intent(in) :: this 
    trunks = this%trunks_i
  end function trunks 

  integer function dead_leaves(this)
    class(fuel_classes_type), intent(in) :: this 
    dead_leaves = this%dead_leaves_i
  end function dead_leaves 

  integer function live_grass(this)
    class(fuel_classes_type), intent(in) :: this 
    live_grass = this%live_grass_i
  end function live_grass 

  logical function moss_classes_present(this)
    class(fuel_classes_type), intent(in) :: this
    moss_classes_present = num_fuel_classes == 8
  end function moss_classes_present

  integer function live_moss(this)
    class(fuel_classes_type), intent(in) :: this
    if (.not. this%moss_classes_present()) then
       write(fates_log(),*) 'live_moss fuel class requested but not present'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    live_moss = this%live_moss_i
  end function live_moss

  integer function dead_moss(this)
    class(fuel_classes_type), intent(in) :: this
    if (.not. this%moss_classes_present()) then
       write(fates_log(),*) 'dead_moss fuel class requested but not present'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    dead_moss = this%dead_moss_i
  end function dead_moss

end module FatesFuelClassesMod
