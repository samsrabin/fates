module FatesEdgeForestMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : num_edge_forest_bins
  use FatesGlobals, only : fates_log
  use FatesGlobals, only : endrun => fates_endrun
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use EDTypesMod, only : ed_site_type
  use FatesPatchMod, only : fates_patch_type
  use FatesEcotypesMod, only : is_patch_forest

  implicit none
  private  ! By default everything is private

  ! Make public necessary subroutines and functions
  public :: calculate_edge_area
  public :: indexx  ! public for unit testing

contains

  ! =====================================================================================

  function get_number_of_forest_patches(site) result(n_forest_patches)
    ! DESCRIPTION
    ! Returns number of forest patches at site
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    !
    ! RETURN VALUE:
    integer :: n_forest_patches
    !
    ! LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch

    n_forest_patches = 0
    currentPatch => site%youngest_patch
    do while(associated(currentPatch))
       if (currentPatch%is_forest) then
          n_forest_patches = n_forest_patches + 1
       end if

       currentPatch => currentPatch%older
    enddo

  end function get_number_of_forest_patches


  function get_number_of_patches(site) result(n_patches)
    ! DESCRIPTION
    ! Returns number of patches at site
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    !
    ! RETURN VALUE:
    integer :: n_patches
    !
    ! LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch

    n_patches = 0
    currentPatch => site%youngest_patch
    do while(associated(currentPatch))
       n_patches = n_patches + 1
       currentPatch => currentPatch%older
    enddo

  end function get_number_of_patches


  subroutine rank_forest_edge_proximity(site, ranks, index_forestpatches_to_allpatches)
    ! DESCRIPTION:
    ! Rank forest patches by their proximity to edge.
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    integer, dimension(:), intent(inout) :: ranks ! Ranks of forest patches (higher = closer to edge)
    integer, dimension(:), intent(inout) :: index_forestpatches_to_allpatches  ! Array with length (number of patches in gridcell), values 0 if not forest and otherwise an index corresponding to which number forest patch this is
    !
    ! LOCAL VARIABLES:
    real(r8), dimension(:), allocatable :: array ! Array to be index-sorted.
    integer :: n_forest_patches  ! Number of patches in above arrays
    type(fates_patch_type), pointer  :: currentPatch
    integer :: f  ! index of current forest patch
    integer :: p  ! index of patch

    ! Skip sites with no forest patches
    n_forest_patches = get_number_of_forest_patches(site)
    if (n_forest_patches == 0) then
       return
    end if

    ! Allocate arrays
    allocate(array(1:n_forest_patches))

    ! Fill arrays
    f = 0
    p = 0
    index_forestpatches_to_allpatches(:) = 0
    currentPatch => site%oldest_patch
    patchloop: do while(associated(currentPatch))
       p = p + 1
       if (.not. currentPatch%is_forest) then
          currentPatch => currentPatch%younger
          cycle
       end if

       f = f + 1
       index_forestpatches_to_allpatches(p) = f

       ! Fill with patch age.
       ! TODO: Add other options. Biomass? Woody biomass?
       array(f) = currentPatch%age

       currentPatch => currentPatch%younger
    end do patchloop

    ! Get indices of sorted forest patches
    call indexx(array, ranks)

    ! Clean up
    deallocate(array)
  end subroutine rank_forest_edge_proximity


! !   subroutine sort_forest_edge_proximity()
! !     ! DESCRIPTION:
! !     ! Sorts a linked list of POINTERS to forest patches based on their proximity to edge.
! !     ! Based on EDCohortDynamicsMod's sort_cohorts().
! !   end subroutine sort_forest_edge_proximity


! !   subroutine insert_forest_patch()
! !     ! DESCRIPTION:
! !     ! Insert forest patch into linked list. Based on EDCohortDynamicsMod's insert_cohort().
! !   end subroutine insert_forest_patch


  subroutine get_fraction_of_forest_in_each_bin(fraction_forest_in_each_bin)
    ! DESCRIPTION:
    ! Get the fraction of forest in each bin.
    ! PLACEHOLDER FOR NOW that just returns 1/num_edge_forest_bins for each.
    ! TODO: Replace this with real equations for each bin as a function of deforested area
    !
    ! ARGUMENTS
    real(r8), dimension(:), pointer, intent(in) :: fraction_forest_in_each_bin
    !
    ! LOCAL VARIABLES
    integer :: b  ! Bin index

    binloop: do b = 1, num_edge_forest_bins
       fraction_forest_in_each_bin(b) = 1._r8 / real(num_edge_forest_bins, r8)
    end do binloop
  end subroutine get_fraction_of_forest_in_each_bin


  subroutine calculate_edge_area(site)
    ! DESCRIPTION:
    ! Loop through forest patches in decreasing order of proximity, calculating the
    ! area of each patch that is in each edge bin.
    !
    ! ARGUMENTS:
    type(ed_site_type), pointer, intent(in) :: site
    !
    ! LOCAL VARIABLES:
    integer, dimension(:), allocatable :: ranks ! Ranks of forest patches (higher = closer to edge)
    integer, dimension(:), allocatable :: index_forestpatches_to_allpatches  ! Array with length (number of patches in gridcell), values 0 if not forest and otherwise an index corresponding to which number forest patch this is
    integer :: n_forest_patches  ! Number of forest patches
    integer :: n_patches  ! Number of patches in site
    real(r8), dimension(num_edge_forest_bins), target :: fraction_forest_in_each_bin

    ! Skip sites with no forest patches
    n_forest_patches = get_number_of_forest_patches(site)
    if (n_forest_patches == 0) then
       return
    end if

    ! Allocate arrays
    allocate(ranks(1:n_forest_patches))
    n_patches = get_number_of_patches(site)
    allocate(index_forestpatches_to_allpatches(1:n_patches))

    ! Get ranks
    call rank_forest_edge_proximity(site, ranks, index_forestpatches_to_allpatches)

    ! Get fraction of forest area in each bin
    call get_fraction_of_forest_in_each_bin(fraction_forest_in_each_bin)

    ! Clean up
    deallocate(ranks)
    deallocate(index_forestpatches_to_allpatches)
  end subroutine calculate_edge_area




  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array.
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx.
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. The conversion to
  ! Fortran 90, and modification to do an index sort was done by Ian Rutt.
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indexx(array, index)

    ! Performs an index sort of \texttt{array} and returns the result in
    ! \texttt{index}. The order of elements in \texttt{array} is unchanged.
    !
    ! This is a GPL-licenced replacement for the Numerical Recipes routine indexx.
    ! It is not derived from any NR code, but are based on a quicksort routine by
    ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    ! in C, and issued under the GNU General Public License. The conversion to
    ! Fortran 90, and modification to do an index sort was done by Ian Rutt.

    real(r8), dimension(:) :: array ! Array to be indexed.
    integer, dimension(:) :: index ! Index of elements of patch_array
    integer :: i

    if (size(array) /= size(index)) then
       write(fates_log(),*) 'ERROR: INDEXX size mismatch.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    else if (size(array) == 0) then
       write(fates_log(),*) 'ERROR: INDEXX array size 0.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do i=1,size(index)
       index(i)=i
    enddo

    call q_sort_index(array,index,1,size(array))

  end subroutine indexx

!==============================================================

  recursive subroutine q_sort_index(numbers,index,left,right)

    !> This is the recursive subroutine actually used by sort_patches.
    !>
    !> This is a GPL-licenced replacement for the Numerical Recipes routine indexx.
    !> It is not derived from any NR code, but are based on a quicksort routine by
    !> Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !> in C, and issued under the GNU General Public License. The conversion to
    !> Fortran 90, and modification to do an index sort was done by Ian Rutt.

    implicit none

    real(r8), dimension(:) :: numbers !> Numbers being sorted
    integer, dimension(:) :: index   !> Returned index
    integer :: left, right           !> Limit of sort region

    integer :: ll,rr
    integer :: pv_int,l_hold, r_hold,pivpos
    real(r8) :: pivot

    ll=left
    rr=right

    l_hold = ll
    r_hold = rr
    pivot = numbers(index(ll))
    pivpos=index(ll)

    do
       if (.not.(ll < rr)) exit

       do
          if  (.not.((numbers(index(rr)) >= pivot) .and. (ll < rr))) exit
          rr=rr-1
       enddo

       if (ll /= rr) then
          index(ll) = index(rr)
          ll=ll+1
       endif

       do
          if (.not.((numbers(index(ll)) <= pivot) .and. (ll < rr))) exit
          ll=ll+1
       enddo

       if (ll /= rr) then
          index(rr) = index(ll)
          rr=rr-1
       endif
    enddo

    index(ll) = pivpos
    pv_int = ll
    ll = l_hold
    rr = r_hold
    if (ll < pv_int)  call q_sort_index(numbers, index,ll, pv_int-1)
    if (rr > pv_int)  call q_sort_index(numbers, index,pv_int+1, rr)

  end subroutine q_sort_index


  end module FatesEdgeForestMod