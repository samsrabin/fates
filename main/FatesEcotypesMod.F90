module FatesEcotypesMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesPatchMod, only : fates_patch_type

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: is_patch_forest

contains

  ! =====================================================================================

  function is_patch_forest(patchptr, forest_tree_fraction_threshold)
  ! DESCRIPTION:
  ! Return boolean: Is this patch "forest"?
  !
  ! ARGUMENTS:
  type(fates_patch_type), intent(in), pointer :: patchptr  ! pointer to patch object
  real(r8), intent(in) :: forest_tree_fraction_threshold ! Tree fraction above which a patch is "forest"
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