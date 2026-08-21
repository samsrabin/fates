#!/usr/bin/env python

#### Standalone generator for the moss (non-vascular phototroph) FATES PFT:
#    adds a 15th FATES PFT column and two new litterclasses (live/dead
#    moss) to build fates_params_moss.json from fates_params_default.json,
#    and adds the fates_vascular parameter to the default file itself.
#    This is deliberately NOT a batch_patch_params.py patch file: that
#    machinery can add a 15th PFT column (via pft_index_swapper.py's numpy
#    indexing) but cannot grow dimensions.fates_litterclass 6 -> 8, nor
#    create a parameter absent from the base file (fates_vascular).
#
#    This script uses only the stdlib json module plus FATES's own
#    write_json helper -- no numpy, no batch_patch_params.py.
#
#    Two independent operations, selected by --mode:
#
#    add-vascular   Add the fates_vascular parameter (all PFTs = 1) to a
#                   FATES parameter JSON that does not yet have it. Used
#                   to add fates_vascular to fates_params_default.json
#                   ahead of the moss PFT work. Safe to add before any
#                   Fortran reads it: FatesTransferParameters claims
#                   parameters by name, so an unclaimed extra parameter is
#                   harmless.
#
#                   fates_params_default.json's on-disk style predates
#                   write_json.py and does not byte-for-byte match what
#                   write_json.traverse_data emits (see _match_repo_style
#                   below for the specifics). Reserializing the whole file
#                   through write_json would therefore reformat every
#                   parameter entry -- a whole-file diff -- which would
#                   defeat the point of this mode being a pure insertion
#                   into an existing, hand-maintained file. So this mode
#                   does NOT reserialize the document: it renders just the
#                   one new parameter block with write_json, fixes up that
#                   fragment's formatting with _match_repo_style, and
#                   splices the result in immediately after the existing
#                   fates_turnover_senleaf_fdrought entry -- the correct
#                   alphabetical slot for fates_vascular, immediately
#                   before fates_wood_density (see sort_parameters.py's
#                   canonical per-dims-group alphabetical ordering). No
#                   existing line in the file is touched.
#
#    make-moss      Build fates_params_moss.json from a base file (which
#                   must already carry fates_vascular, i.e. the output of
#                   --mode add-vascular). This is a brand-new file, so it
#                   is written in full through write_json.traverse_data
#                   and then passed through the same _match_repo_style
#                   fixup, so its on-disk formatting matches every other
#                   FATES parameter file (write_json's raw output does
#                   not -- see _match_repo_style).
#
#                   attributes.history in the output records the actual
#                   invocation (sys.argv, no timestamp -- see
#                   _moss_history_text), so the committed file's bytes
#                   depend on exactly how this script was called: path
#                   spellings, the --fout name, argument order. To
#                   reproduce parameter_files/fates_params_moss.json as
#                   committed byte-for-byte, run this canonical command
#                   from the src/fates directory, verbatim, with the
#                   relative paths shown (an absolute path here would
#                   leak this machine's directory layout into a
#                   committed FATES parameter file):
#
#                       tools/make_moss_params.py --mode make-moss \
#                           --fin parameter_files/fates_params_default.json \
#                           --fout parameter_files/fates_params_moss.json

import argparse
import io
import json
import os
import re
import sys

# 1. Get the directory of the currently executing script
script_dir = os.path.dirname(os.path.abspath(__file__))
# 2. Insert this directory at the beginning of the Python search path so
#    "import write_json" finds the sibling module regardless of cwd.
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

import write_json


# ----------------------------------------------------------------------
# On-disk style fixup
# ----------------------------------------------------------------------

# Matches a nested-dict's opening line as emitted by write_json.traverse_data
# (e.g. '  "attributes":{'), which has no space before the brace.
_TIGHT_BRACE_RE = re.compile(r'":\{$', re.M)

# Matches the top-level "attributes" dict as rendered (multi-line) by
# write_json.traverse_data, so it can be collapsed onto one line to match
# every other FATES parameter file's layout. Attribute values are JSON
# strings, so any literal newline inside one is escaped ("\n"), not raw --
# the only raw newlines in this span are the structural ones this regex
# is written to find.
_ATTRIBUTES_RE = re.compile(r'  "attributes": \{\n(.*?)\n  \},\n', re.S)


def _match_repo_style(text):
    """Reformat write_json.traverse_data's raw output to match the on-disk
    conventions used by every other FATES parameter file (e.g.
    fates_params_default.json):

    - a space before each nested-dict's opening brace ('"key": {', not
      write_json's raw '"key":{');
    - no space after the comma between a 2-D "data" array's outer
      sublists ('],[', not write_json's/json.dumps's raw '], [');
    - the top-level "attributes" dict collapsed onto one line (write_json
      always expands a dict multi-line; every FATES parameter file keeps
      "attributes" on a single line).

    Used by both --mode add-vascular (applied to its one new-parameter
    fragment; the "attributes" substitution is simply a no-op there,
    since the fragment never contains an "attributes" key) and --mode
    make-moss (applied to the whole rendered document).
    """
    text = _TIGHT_BRACE_RE.sub('": {', text)
    text = text.replace("], [", "],[")

    def _collapse_attributes(match):
        items = [line.strip().rstrip(",") for line in match.group(1).split("\n")]
        return '  "attributes": {' + ", ".join(items) + "},\n"

    text = _ATTRIBUTES_RE.sub(_collapse_attributes, text)
    return text


# ----------------------------------------------------------------------
# Shared constants
# ----------------------------------------------------------------------

GRASS_PFT_INDEX = 11  # 0-based index of arctic_c3_grass (FATES PFT 12)
MOSS_PFT_NAME = "non_vascular_phototroph"  # NVP's name for the moss PFT

# Values applied to the moss (new, 15th) FATES PFT column after it is
# seeded by copying the arctic_c3_grass column, harvested from NVP's
# moss column (8382939b9:parameter_files/fates_params_default_moss.json)
# unless noted otherwise. See the "harvested" and "corrections" section
# markers below for which values come straight from that file versus
# are deliberately overridden here.
#
# NOTE: fates_allom_fnrt_prof_mode is intentionally ABSENT from this
# dict. NVP's moss column sets it to 4, a no-roots profile that exists
# only on that branch; our FATES supports modes 1-3, and moss
# transpiration has to be extracted through a real root profile or the
# water budget does not close -- see fates_allom_fnrt_prof_a below for
# how that profile is shaped. Moss keeps the grass-copied mode 3
# instead.
MOSS_PFT_OVERRIDES = {
    # --- harvested from NVP's moss column: taken as-is because staying
    #     aligned with that branch is the point ---
    # NVP's name for the moss PFT, kept verbatim.
    "fates_pftname": MOSS_PFT_NAME,
    # Porada et al. (2013)-based moss photosynthetic capacity (grass is
    # 86.0); harvested from NVP.
    "fates_leaf_vcmax25top": 30.0,  # dims include fates_leafage_class too
    # NVP's value, harvested as-is; no opinion offered about it here.
    "fates_rad_leaf_clumping_index": 10.0,
    # fates_leaf_slatop and fates_woody (below) are no-ops: NVP's harvested
    # moss value for each already equals the arctic_c3_grass value the
    # moss column is seeded from. Listed anyway for completeness/parity
    # with NVP's column -- not a moss-vs-grass difference (slatop's
    # provenance is still Porada et al. 2013-based SLA, same source as
    # vcmax25top above). fates_woody = 0 is also what routes moss
    # through the non-woody paths: live biomass to the live fuel pool,
    # stem litter to leaf fines, no treefall disturbance, dbh forced
    # from leaf carbon.
    "fates_leaf_slatop": 0.027,
    "fates_woody": 0,
    # Moss has no stomata; the moss CO2 path replaces the stomatal solve
    # with boundary-layer diffusion, so these three are unused for moss.
    # Zeroed so any stray use shows up as zero rather than as a
    # plausible number.
    "fates_leaf_stomatal_intercept": 0.0,
    "fates_leaf_stomatal_slope_ballberry": 0.0,
    "fates_leaf_stomatal_slope_medlyn": 0.0,
    # btran_on_ag_none (biogeophys/LeafBiophysicsMod.F90:176): do not
    # apply btran to vcmax or jmax. Moss capacity is limited instead by
    # the fwet proxy (vcmax * min(1, fwet/0.6)); applying btran on top
    # would double-count water limitation.
    "fates_leaf_agross_btran_model": 0,
    # ievergreen (main/FatesConstantsMod.F90:75; grass is 2). Moss
    # carries its thallus year-round, and evergreen also sidesteps the
    # deciduous phenflush_fraction requirement enforced at
    # main/EDPftvarcon.F90:1219.
    "fates_phen_leaf_habit": 1,
    # 0.01 is near-opaque; harvested from NVP.
    "fates_rad_leaf_taunir": 0.01,
    "fates_rad_leaf_tauvis": 0.01,
    "fates_rad_stem_taunir": 0.01,
    "fates_rad_stem_tauvis": 0.01,
    # Leaf-angle orientation index; 0 is random/spherical (grass is
    # -0.23). Harvested from NVP; no rationale is recorded on that side
    # for why moss should have it this way.
    "fates_rad_leaf_xl": 0.0,
    # The flag that identifies moss -- grass and every other PFT are 1.
    # Nothing reads it yet; the Fortran that does arrives in the next
    # task.
    "fates_vascular": 0,
    # --- corrections applied here: for these three, NVP's moss column
    #     still holds the grass values, and we deliberately override
    #     them ---
    # Reproductive allocation is two branches
    # (parteh/PRTAllometricCarbonMod.F90:1074-1078): below
    # dbh_repro_threshold, repro_fraction = seed_alloc; above it,
    # repro_fraction = seed_alloc + seed_alloc_mature. Grass survives
    # seed_alloc = 0.0 because it clears its 3.0 cm threshold and then
    # collects seed_alloc_mature = 0.25. Moss (~0.03 cm dbh) would
    # never clear a 3.0 cm threshold, so it would sit on the immature
    # branch at 0.0 forever: permanent zero reproductive allocation.
    # Dropping the threshold to 0.001 is therefore the whole fix -- it
    # puts moss on the mature branch, where the inherited 0.25
    # applies. Both seed_alloc (0.0) and seed_alloc_mature (0.25) are
    # inherited from the grass copy and deliberately NOT overridden --
    # removing seed_alloc_mature later would silently sterilize moss.
    "fates_recruit_seed_dbh_repro_threshold": 0.001,
    # fates_allom_dbh_maxheight is inherited from grass (20.0),
    # deliberately NOT overridden: the TRS tree test
    # (allom_dbh_maxheight > min_max_dbh_for_trees, 15 cm;
    # main/FatesConstantsMod.F90:134) is irrelevant here because every
    # TRS gate is conjoined with an hlm_regeneration_model test
    # (parteh/PRTAllometricCarbonMod.F90:1069-1070,1083-1084), and
    # fates_regeneration_model defaults to default -- this project does
    # not run the experimental TRS. (If TRS were ever switched on,
    # moss would be classified a tree.) The trade-off: this parameter
    # is also the diameter at which height and max leaf biomass
    # saturate (d2h_2pwr's h = p1*min(d,dbh_maxh)**p2 and
    # dh2blmax_3pwr_grass's duse = min(d,dbh_maxh);
    # biogeochem/FatesAllometryMod.F90), so inheriting 20.0 leaves moss
    # height effectively unbounded under grass's power-law allometry --
    # it saturates only around 1.23 m, where a moss-specific 0.1 cm
    # would have capped it near 4.2 cm (moss recruits at 2 cm height,
    # about 0.032 cm dbh). The mat-thickness height allometry mode is
    # the principled fix; this is an accepted limitation recorded in
    # the design spec.
    # A realistic 2 cm moss recruit -- bounds only the recruit's
    # initial size, not the (now grass-inherited, unbounded) maximum
    # size set by fates_allom_dbh_maxheight above. Grass's 0.11 m would
    # inflate the allometric per-plant target biomass. Raise it only if
    # cohort termination floors actually cull moss in testing (watch
    # the FATES_MORTALITY_TERMINATION_* history variables).
    "fates_recruit_height_min": 0.02,
    # Concentrates the rooting profile in soil layer 1 so moss water
    # status tracks surface moisture. Needed because we keep mode 3
    # rather than NVP's no-roots mode 4 (see the NOTE above); moss
    # transpiration is extracted through this profile, and an all-zero
    # one would break the water budget.
    "fates_allom_fnrt_prof_a": 30.0,
    # fates_allom_fnrt_prof_mode is deliberately NOT overridden: it keeps
    # the grass-copied value of 3 (see the NOTE above).
}

DEAD_LEAVES_INDEX = 4  # 0-based index of "dead leaves" in fates_litterclass
# All ten litterclass-dimensioned parameters grow by copying the
# "dead leaves" values into both new (7th, 8th) slots, except:
LITTERCLASS_OVERRIDE = {
    "fates_frag_maxdecomp": [999.0, 1.0],  # index 6 unread; index 7 read
    "fates_litterclass_name": ["live moss", "dead moss"],
}

MOSS_ATTRIBUTE_KEY = "moss_pft_caveat"
MOSS_ATTRIBUTE_TEXT = (
    "fates_hlm_pft_map row 4 (HLM natpft index 4) is remapped entirely to "
    "FATES PFT 15 (non_vascular_phototroph) in this file. This parameter "
    "file is therefore only sensible for surface datasets where HLM PFT 4 "
    "(broadleaf_evergreen_tropical_tree) carries no real area -- e.g. at "
    "a tropical site it would silently turn broadleaf evergreen tropical "
    "tree area into moss."
)

# Replaces the inherited attributes.history line (which reads "...copied
# from: ../parameter_files/fates_params_default.cdl." -- untrue of this
# generated file) with the actual invocation that produced this file,
# following the change_log convention in sort_parameters.py:97-102 and
# batch_patch_params.py:30 (both build ' '.join(sys.argv[1:])), with two
# deliberate deviations: no datetime.now() prefix -- a timestamp would
# make every regeneration produce different bytes, defeating the
# byte-reproducibility this project verifies every round -- and the
# whole command (sys.argv[0] plus every argument, not just
# sys.argv[1:]) is kept, so the invocation is reproducible from the
# string alone.
def _moss_history_text():
    return " ".join(sys.argv)


# ----------------------------------------------------------------------
# Helpers for walking the innermost (fates_pft) axis of nested lists
# ----------------------------------------------------------------------

def _append_last_axis(data, source_index):
    """Recursively append data[...][source_index] as a new last element
    along the innermost axis of a (possibly nested) list, in place."""
    if isinstance(data, list) and data and isinstance(data[0], list):
        for row in data:
            _append_last_axis(row, source_index)
    elif isinstance(data, list):
        data.append(data[source_index])
    else:
        raise TypeError(f"Expected a list, got {type(data)}")


def _set_last_element(data, value):
    """Recursively set the last element along the innermost axis of a
    (possibly nested) list to value, in place."""
    if isinstance(data, list) and data and isinstance(data[0], list):
        for row in data:
            _set_last_element(row, value)
    elif isinstance(data, list):
        data[-1] = value
    else:
        raise TypeError(f"Expected a list, got {type(data)}")


# ----------------------------------------------------------------------
# --mode add-vascular
# ----------------------------------------------------------------------

def _fates_vascular_entry(npft):
    return {
        "dtype": "integer",
        "dims": ["fates_pft"],
        "long_name": "Whether plant is vascular",
        "units": "unitless - 1 true, 0 false",
        "data": [1] * npft,
    }


def do_add_vascular(args):
    with open(args.inputfname, "r") as fin:
        raw_text = fin.read()
    data = json.loads(raw_text)

    if "fates_vascular" in data["parameters"]:
        existing = data["parameters"]["fates_vascular"]["data"]
        if all(x == 1 for x in existing):
            # Already added (e.g. a repeat run); write the input back out
            # unchanged so this mode is idempotent / reproducible.
            with open(args.outputfname, "w") as fout:
                fout.write(raw_text)
            if not args.silent:
                print("fates_vascular already present (all 1); no change made.")
            return
        raise SystemExit(
            "fates_vascular is already present but is not all 1; refusing "
            "to overwrite blindly -- check the input file by hand."
        )

    npft = data["dimensions"]["fates_pft"]
    entry = _fates_vascular_entry(npft)

    # Render just the new parameter block with FATES's own writer, so its
    # internal formatting (the dtype/dims/long_name/units/data lines)
    # matches every other entry in the file byte-for-byte. last_att=False
    # because it will be inserted ahead of an existing entry (it is not
    # the last parameter in the file), so it needs a trailing comma.
    buf = io.StringIO()
    write_json.traverse_data(
        buf, entry, indent_level=2, current_key="fates_vascular", last_att=False
    )
    fragment = _match_repo_style(buf.getvalue())

    # Splice the fragment in immediately after the existing
    # fates_turnover_senleaf_fdrought entry -- the correct alphabetical
    # slot for fates_vascular (immediately before fates_wood_density; see
    # sort_parameters.py's canonical per-dims-group alphabetical order) --
    # so the diff is a pure insertion -- no existing line is modified.
    # (fates_vascular's metadata is still modeled on fates_woody's, above;
    # only the splice location is fates_turnover_senleaf_fdrought.)
    anchor_open = '    "fates_turnover_senleaf_fdrought": {\n'
    if raw_text.count(anchor_open) != 1:
        raise SystemExit(
            f"Expected exactly one '{anchor_open.strip()}' line; found "
            f"{raw_text.count(anchor_open)}. Refusing to splice blindly."
        )
    idx_open = raw_text.index(anchor_open)
    close_line = "    },\n"
    idx_close = raw_text.index(close_line, idx_open)
    insertion_point = idx_close + len(close_line)

    new_text = raw_text[:insertion_point] + fragment + raw_text[insertion_point:]

    # Sanity check: the spliced text must be valid JSON with the expected
    # shape before we write it out.
    check = json.loads(new_text)
    if check["parameters"]["fates_vascular"]["data"] != [1] * npft:
        raise SystemExit("Post-splice sanity check failed unexpectedly.")

    with open(args.outputfname, "w") as fout:
        fout.write(new_text)
    if not args.silent:
        print(f"Added fates_vascular (all 1, length {npft}) to {args.outputfname}")


# ----------------------------------------------------------------------
# --mode make-moss
# ----------------------------------------------------------------------

def do_make_moss(args):
    with open(args.inputfname, "r") as fin:
        data = json.load(fin)

    params = data["parameters"]
    if "fates_vascular" not in params:
        raise SystemExit(
            "Base file has no fates_vascular parameter; this mode only "
            "overrides parameters that already exist on the base file, so "
            "run --mode add-vascular on it first."
        )

    npft_old = data["dimensions"]["fates_pft"]
    if npft_old != 14:
        raise SystemExit(f"Expected 14 base FATES PFTs, found {npft_old}")
    nlitt_old = data["dimensions"]["fates_litterclass"]
    if nlitt_old != 6:
        raise SystemExit(f"Expected 6 base litterclasses, found {nlitt_old}")
    if params["fates_litterclass_name"]["data"][DEAD_LEAVES_INDEX] != "dead leaves":
        raise SystemExit(
            "fates_litterclass_name[4] is not 'dead leaves'; litterclasses "
            "were reordered"
        )

    # --- 1. Append the 15th (moss) FATES PFT column, seeded from the
    #        arctic_c3_grass column, for every fates_pft-dimensioned
    #        parameter except fates_hlm_pft_map (handled separately: an
    #        area-fraction map is not a "starting point" template). -----
    for key, entry in params.items():
        if "fates_pft" not in entry["dims"]:
            continue
        if entry["dims"][-1] != "fates_pft":
            raise SystemExit(f"{key}: expected fates_pft as the last dim")
        if key == "fates_hlm_pft_map":
            continue
        _append_last_axis(entry["data"], GRASS_PFT_INDEX)

    data["dimensions"]["fates_pft"] = npft_old + 1

    # --- 2. Apply the moss-specific PFT parameter overrides -------------
    missing = sorted(set(MOSS_PFT_OVERRIDES) - set(params))
    if missing:
        raise SystemExit(f"Base file is missing override targets: {missing}")
    for key, value in MOSS_PFT_OVERRIDES.items():
        _set_last_element(params[key]["data"], value)

    # --- 3. fates_hlm_pft_map: the new column starts at 0 for every row,
    #        then HLM PFT 4 (row index 3) is fully remapped to moss.
    #        Row 4 is rebuilt unconditionally, not patched: HLM PFT 4 is
    #        assigned to moss by this file's design, so the rebuild
    #        itself does not depend on the base file's prior content.
    #        But before rebuilding, the base file's row 4 is validated
    #        against that same assumption (all of HLM PFT 4 -> FATES PFT
    #        1) -- a surprising row 4 means the assumption this script is
    #        built on may no longer hold, which is a hard error, not
    #        something to silently paper over. The row-sum loop below
    #        then validates every row, including the just-rebuilt row 4.
    hlm_map = params["fates_hlm_pft_map"]["data"]
    for row in hlm_map:
        row.append(0.0)

    # Hard error, not a note: if the base file's row 4 isn't (within the
    # same 1e-12 tolerance as the row-sum check below) "all of HLM PFT 4
    # -> FATES PFT 1", the assumption this generator is built on may no
    # longer hold. That is a deliberate decision for a person to make --
    # revisit the moss remap target in this script -- not something to
    # force through.
    old_row4 = hlm_map[3]
    expected_row4 = [1.0] + [0.0] * (len(old_row4) - 1)
    if any(abs(o - e) > 1e-12 for o, e in zip(old_row4, expected_row4)):
        raise SystemExit(
            f"fates_hlm_pft_map row 4 is {old_row4}, not entirely FATES "
            "PFT 1 (expected [1.0, 0.0, ..., 0.0]). This generator "
            "requires HLM PFT 4 to map entirely to FATES PFT 1, because it "
            "reassigns that slot to moss. If the base file's row 4 has "
            "legitimately changed upstream, the moss remap target in this "
            "script needs to be revisited -- do not force this check to "
            "pass."
        )

    hlm_map[3] = [0.0] * (npft_old + 1)
    hlm_map[3][npft_old] = 1.0
    for irow, row in enumerate(hlm_map):
        if abs(sum(row) - 1.0) > 1e-12:
            raise SystemExit(f"hlm_pft_map row {irow+1} sums to {sum(row)}, not 1.0")

    # --- 4. Grow fates_litterclass 6 -> 8: slots 7-8 (live/dead moss)
    #        copy the "dead leaves" values as starting points, except
    #        fates_frag_maxdecomp and fates_litterclass_name ------------
    for key, entry in params.items():
        if "fates_litterclass" not in entry["dims"]:
            continue
        if entry["dims"] != ["fates_litterclass"]:
            raise SystemExit(f"{key}: unexpected fates_litterclass dims")
        if key in LITTERCLASS_OVERRIDE:
            entry["data"].extend(LITTERCLASS_OVERRIDE[key])
        else:
            dead_leaves_value = entry["data"][DEAD_LEAVES_INDEX]
            entry["data"].extend([dead_leaves_value, dead_leaves_value])

    data["dimensions"]["fates_litterclass"] = nlitt_old + 2

    # --- 5. Record the HLM-PFT-4 caveat, and replace the (now-false)
    #        inherited history line with accurate, static provenance -----
    data["attributes"][MOSS_ATTRIBUTE_KEY] = MOSS_ATTRIBUTE_TEXT
    data["attributes"]["history"] = _moss_history_text()

    # This is a brand-new file (no prior on-disk formatting to preserve),
    # so it is written in full through write_json, then passed through
    # _match_repo_style so its on-disk formatting matches every other
    # FATES parameter file (write_json's raw output does not).
    buf = io.StringIO()
    write_json.traverse_data(buf, data)
    text = _match_repo_style(buf.getvalue())
    with open(args.outputfname, "w") as fout:
        fout.write(text)
    if not args.silent:
        print(f"Wrote {args.outputfname}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate FATES moss (non_vascular_phototroph) parameter "
            "data: add fates_vascular to a base parameter file "
            "(--mode add-vascular), or build the moss parameter file from "
            "a base file that already has it (--mode make-moss). See "
            "this script's module docstring for details."
        )
    )
    parser.add_argument(
        "--mode",
        dest="mode",
        choices=["add-vascular", "make-moss"],
        required=True,
        help=(
            "add-vascular: add the fates_vascular parameter (all 1) to a "
            "base file that lacks it. make-moss: build the moss parameter "
            "file from a base file that already has fates_vascular."
        ),
    )
    parser.add_argument(
        "--fin", "--input", dest="inputfname", type=str,
        help="Input filename. Required.", required=True,
    )
    parser.add_argument(
        "--fout", "--output", dest="outputfname", type=str,
        help="Output filename. Required.", required=True,
    )
    parser.add_argument(
        "--silent", dest="silent", help="Suppress output messages",
        action="store_true",
    )

    args = parser.parse_args()

    if args.mode == "add-vascular":
        do_add_vascular(args)
    else:
        do_make_moss(args)


if __name__ == "__main__":
    main()
