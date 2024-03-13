module FatesEcotypesMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesPatchMod, only : fates_patch_type

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: is_patch_forest
  
  real(r8) :: forest_tree_fraction_threshold = 0.10_r8  ! FAO: "Forest" is > 10% tree cover

contains

  ! =====================================================================================

  function is_patch_forest(patchptr)
  ! DESCRIPTION:
  ! Return boolean: Is this patch "forest"?
  !
  ! ARGUMENTS:
  type(fates_patch_type), intent(in), pointer :: patchptr  ! pointer to patch object
  !
  ! RETURN VALUE
  logical :: is_patch_forest
  !
  ! LOCAL VARIABLES
  real(r8) :: tree_fraction = 0._r8

  if (patchptr%area > 0._r8) then
      tree_fraction = patchptr%total_tree_area / patchptr%area
  end if

  is_patch_forest = tree_fraction > forest_tree_fraction_threshold

  end function is_patch_forest


end module FatesEcotypesMod