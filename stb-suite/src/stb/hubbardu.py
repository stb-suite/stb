#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

# 2.2.0: numbered-section report (matching the other Stage 1 prep tools);
# -s/--structure and -c/--calc now default to 'structure.fdf'/'calc.fdf';
# --save-report; the DFT+U perturbation block (and a new, previously-missing
# single-point-SCF enforcement -- the reference run was never forced off a
# relaxation-style --calc template, silently invalidating the fixed-geometry
# assumption the whole linear-response method depends on) now go into a
# config_extra.fdf sidecar + '%include', the same pattern strain.py/
# elastic_inputs.py/phonons_create.py already use, instead of being embedded
# directly inline in calc.fdf.
VERSION = "2.2.0"

import os
import re
import sys
import json
import shutil
import argparse
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.dftu_data import SHELL_NAMES, DEFAULT_SHELL, ldau_proj_block as _base_ldau_proj_block
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source
from stb.core.symmetry import find_inequivalent_sites

_LABEL_RE = re.compile(r'^\s*SystemLabel\s+(\S+)', re.IGNORECASE | re.MULTILINE)
_EXISTING_LDAU_RE = re.compile(r'LDAU\.proj|LDAU\.PotentialShift|DFTU\.proj|DFTU\.PotentialShift', re.IGNORECASE)

REPORT_FILE = "hubbardu_stage1.txt"
EXTRA_FDF_FILE = "config_extra.fdf"
TEMPLATE_STRUCTURE_NAME = "_template_structure.fdf"
TEMPLATE_CALC_NAME = "_template_calc.fdf"
PSEUDO_SNAPSHOT_DIRNAME = "_pseudopotentials"

# Forces every run folder (reference here, and scf_alpha_*/frozen_alpha_* in
# stb-hubbarduAlphas, which also calls write_run_folder() below) to a pure
# single-point SCF, regardless of what --calc itself sets -- the Cococcioni-
# de Gironcoli linear-response method requires every run at the EXACT SAME
# geometry; letting a relaxation-style --calc template move the ions would
# silently invalidate the comparison between the reference and every
# perturbed alpha.
SINGLE_POINT_BLOCK = (
    "# Forces a pure single-point SCF (no ionic or cell relaxation) -- the\n"
    "# Cococcioni-de Gironcoli linear-response method requires every run\n"
    "# (reference and every perturbed alpha) at the EXACT SAME geometry;\n"
    "# letting --calc's own relaxation settings move the ions would silently\n"
    "# invalidate the comparison, regardless of what --calc itself sets.\n"
    "MD.TypeOfRun       CG\n"
    "MD.Steps           0\n"
    "MD.NumCGsteps      0\n"
    "MD.VariableCell    false\n"
)


def get_system_label(calc_text, structure_text=None):
    """SystemLabel, searched in calc_text first (the common case: it's a
    plain fdf keyword, not something SIESTA allows only in an %include'd
    file) and, if not found there, in structure_text too -- some templates
    put every fdf keyword in structure.fdf and just %include it verbatim.
    Falls back to "siesta" (SIESTA's own default) only if neither has it.
    """
    match = _LABEL_RE.search(calc_text)
    if match:
        return match.group(1)
    if structure_text:
        match = _LABEL_RE.search(structure_text)
        if match:
            return match.group(1)
    return "siesta"


def ldau_perturbation_block(species, n, l, u, j, potential_shift):
    """The %block LDAU.proj stanza plus, when potential_shift, the
    LDAU.PotentialShift flag that makes SIESTA treat `u` as a rigid
    perturbation on the shell (Cococcioni & de Gironcoli, PRB 71, 035105,
    2005) instead of the real Hubbard term -- and, as a side effect, makes it
    print the shell's occupation to the .out file (verified against SIESTA's
    own source, Src/dftu.F, where that print statement is gated on this exact
    flag).

    Named differently from core.dftu_data.ldau_proj_block (which this wraps)
    on purpose: that one takes a list of {species,n,l,u,j} dicts and never
    adds PotentialShift, so the two are NOT interchangeable despite doing
    related things -- a same-named pair here previously made it easy to
    import the wrong one and silently lose the PotentialShift line.
    """
    block = _base_ldau_proj_block([{"species": species, "n": n, "l": l, "u": u, "j": j}])
    if potential_shift:
        return "LDAU.PotentialShift T\n" + block
    return block


def write_run_folder(output_dir, folder_name, calc_template, structure_path, extra_fdf_text):
    """Writes structure.fdf + calc.fdf + config_extra.fdf into a new run
    folder. `extra_fdf_text` (the DFT+U perturbation block, plus whatever
    else a caller adds -- e.g. stb-hubbarduAlphas' own DM.UseSaveDM/
    MaxSCFIterations lines) is combined with SINGLE_POINT_BLOCK above and
    written to config_extra.fdf; calc.fdf itself is then just '%include
    config_extra.fdf' followed by the untouched original template
    (core.structure_io.prepend_include) -- NOT the override text embedded
    directly inline, matching strain.py/elastic_inputs.py/phonons_create.py's
    own convention (extracted here once this became a case of the same
    "SIESTA's fdf reader is first-occurrence-wins, so an override must
    physically precede any pre-existing occurrence of the same directive in
    the base template" need those tools already solved this same way).

    If the folder already exists (e.g. re-running stb-hubbarduAlphas with
    different --alphas/--frozen-iterations on the same --dir), any stale
    calc.out from a PREVIOUS calc.fdf is removed -- it belongs to input that
    no longer exists on disk, and leaving it would let stb-hubbarduAnalysis
    silently parse a result that doesn't match the regenerated input.

    Returns (folder, combined_extra_text) -- the text is returned so a
    caller's own report can preview exactly what was written (mirrors
    phonons_create.py's own config_extra.fdf preview).
    """
    folder = os.path.join(output_dir, folder_name)
    os.makedirs(folder, exist_ok=True)
    stale_output = os.path.join(folder, "calc.out")
    if os.path.exists(stale_output):
        os.remove(stale_output)
    shutil.copy(structure_path, os.path.join(folder, "structure.fdf"))
    combined_extra = (
        "# Auto-generated by stb-hubbardu -- DFT+U perturbation "
        "(Cococcioni-de Gironcoli) + single-point SCF enforcement.\n"
        + extra_fdf_text + "\n" + SINGLE_POINT_BLOCK)
    with open(os.path.join(folder, EXTRA_FDF_FILE), 'w') as f:
        f.write(combined_extra)
    with open(os.path.join(folder, "calc.fdf"), 'w') as f:
        f.write(structure_io.prepend_include(calc_template, EXTRA_FDF_FILE))
    return folder, combined_extra


def check_no_existing_ldau(calc_template, calc_path):
    """Errors out (rather than silently doubling up the block) if the user's
    own template already sets up DFT+U -- e.g. re-running stb-hubbardu by
    mistake on a template that was itself generated by a previous run of this
    same workflow. Checks both LDAU.* and the DFTU.* spelling SIESTA accepts
    as a synonym for the same keywords. Also called from stb-hubbarduAlphas
    against the saved template snapshot, in case a user hand-edits it (e.g.
    to fix an SCF convergence problem) and reintroduces one of these blocks.

    Since write_run_folder() now puts the DFT+U block into a config_extra.fdf
    sidecar (not embedded inline in calc.fdf), the LDAU/DFTU regex alone can
    no longer catch the single most likely real mistake this check exists
    for: pointing --calc straight at a PREVIOUSLY-GENERATED run folder's own
    calc.fdf (e.g. 'hubbardu_runs/reference/calc.fdf'), which today is just a
    plain '%include config_extra.fdf' line -- no LDAU/DFTU keyword in sight,
    even though it very much already has DFT+U set up one file over. A
    genuine, never-processed template a user hands-writes should essentially
    never already reference this suite's own auto-generated sidecar
    filename, so that '%include' line alone is treated as equally
    disqualifying, regardless of which prep tool (stb-hubbardu/
    stb-hubbarduAlphas/stb-elasticInputs/stb-strain/stb-phononsCreate all use
    the same EXTRA_FDF_FILE convention) actually wrote it.
    """
    if _EXISTING_LDAU_RE.search(calc_template):
        print(color_text(
            f"Error: '{calc_path}' already contains an LDAU.proj/DFTU.proj block or "
            "LDAU.PotentialShift/DFTU.PotentialShift -- this looks like it's already set up "
            "for DFT+U (possibly a previous stb-hubbardu output). Point --calc at "
            "a plain template without any DFT+U setup.", 'red'))
        sys.exit(1)
    if f"%include {EXTRA_FDF_FILE}" in calc_template:
        print(color_text(
            f"Error: '{calc_path}' already '%include {EXTRA_FDF_FILE}' -- a genuine, "
            "unprocessed template should never reference this suite's own auto-generated "
            "sidecar file. This looks like the calc.fdf a previous stb-hubbardu/"
            "stb-hubbarduAlphas run already wrote (its DFT+U perturbation now lives in a "
            f"sibling '{EXTRA_FDF_FILE}', not inline). Point --calc at your own original "
            "template instead.", 'red'))
        sys.exit(1)


# Same set recognized elsewhere in the suite (translate.py's readers, clean.py's
# --keep default, cohesive_energy.py/phonons_create.py's pseudo lookups).
PSEUDO_EXTENSIONS = (".psf", ".psml")

# Files this workflow writes/manages itself in every run folder -- never let a
# --pseudo-dir that happens to also contain leftovers from a previous SIESTA
# run (calc.fdf, structure.fdf, siesta.DM, siesta.out, CLOCK, MESSAGES, ...)
# silently overwrite the just-generated, alpha-specific versions. Filtering by
# PSEUDO_EXTENSIONS above already prevents this in practice; this is a second,
# independent guard so a future extension added to PSEUDO_EXTENSIONS (or a
# same-named pseudopotential-looking file) still can't clobber them.
MANAGED_FILENAMES = frozenset({"calc.fdf", "structure.fdf", "run_manifest.json",
                                EXTRA_FDF_FILE, TEMPLATE_STRUCTURE_NAME, TEMPLATE_CALC_NAME})


def copy_pseudopotentials(src_dir, dest_dir, species=None):
    """Copies only recognized pseudopotential files (PSEUDO_EXTENSIONS,
    case-insensitive; non-recursive) from src_dir into dest_dir. Silently
    skips anything else -- src_dir is often the user's working directory
    rather than a pseudopotentials-only folder, and it must never be allowed
    to drag along stray copies of calc.fdf/structure.fdf/a .DM/SIESTA output
    files that would overwrite the ones this workflow just generated (see
    MANAGED_FILENAMES). Returns the list of filenames copied.

    If `species` is given (an iterable of chemical symbols), only files
    whose base name (before any '-'/'_' suffix, e.g. the "Ga" in
    "Ga-gga.nosemicore.psf") case-insensitively matches one of them are
    copied. src_dir is frequently one of the bundled banks (BANKS,
    ~70-80 elements each) -- without this filter every element in the bank
    gets copied into the run folder instead of just the ones the structure
    actually needs.
    """
    wanted = {s.lower() for s in species} if species is not None else None
    copied = []
    for name in sorted(os.listdir(src_dir)):
        if name in MANAGED_FILENAMES or name.endswith(".DM"):
            continue
        if not name.lower().endswith(PSEUDO_EXTENSIONS):
            continue
        if wanted is not None:
            stem = name
            for ext in PSEUDO_EXTENSIONS:
                if stem.lower().endswith(ext):
                    stem = stem[:-len(ext)]
                    break
            symbol = re.split(r'[-_]', stem, maxsplit=1)[0]
            if symbol.lower() not in wanted:
                continue
        src = os.path.join(src_dir, name)
        if os.path.isfile(src):
            shutil.copy(src, os.path.join(dest_dir, name))
            copied.append(name)
    return copied


def alias_pseudopotential_files(dest_dir, species, alias_label):
    """Duplicates every already-copied pseudopotential file for `species` in
    `dest_dir` under `alias_label` instead (e.g. 'Mn.psml' -> 'Mn_pert.psml',
    'Mn-gga.nosemicore.psf' -> 'Mn_pert-gga.nosemicore.psf') -- SIESTA looks
    up a species' pseudopotential by its exact species label, so the aliased
    species (same Z, different label -- see structure_io.alias_single_atom_
    species) needs its own copy under the new name, not just the original
    element's. Returns the list of new filenames written.
    """
    written = []
    for name in sorted(os.listdir(dest_dir)):
        if not name.lower().endswith(PSEUDO_EXTENSIONS):
            continue
        ext = next(e for e in PSEUDO_EXTENSIONS if name.lower().endswith(e))
        stem = name[: len(name) - len(ext)]
        symbol = re.split(r'[-_]', stem, maxsplit=1)[0]
        if symbol.lower() != species.lower():
            continue
        suffix = stem[len(symbol):]
        new_name = f"{alias_label}{suffix}{name[len(stem):]}"
        shutil.copy(os.path.join(dest_dir, name), os.path.join(dest_dir, new_name))
        written.append(new_name)
    return written


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 3: generates the 'reference' folder for a Hubbard U linear-response calculation (Cococcioni & de Gironcoli, Phys. Rev. B 71, 035105, 2005).", 'bold')}
Workflow: (1) this tool generates 'reference/' -- run SIESTA in it yourself;
(2) once it's converged, run stb-hubbarduAlphas to generate the perturbed
'scf_alpha_*'/'frozen_alpha_*' folders (it copies the reference .DM in for
you); (3) run SIESTA in those, then stb-hubbarduAnalysis to compute U.

If 'reference/' doesn't converge: DFT+U SCF on correlated oxides is often
slow/oscillatory with default settings. Things worth trying in your own
calc.fdf template: a lower DM.MixingWeight (~0.05), a higher DM.NumberPulay
(~6), and a looser SCF.H.Tolerance (~0.01 eV).

If --species appears more than once in --structure, --atom-index is required:
SIESTA's %%block LDAU.proj perturbs an entire species label at once, so with
no --atom-index every atom of that species would be perturbed simultaneously
-- measuring a collective response instead of the single isolated-atom
response the Cococcioni-de Gironcoli supercell method requires. When needed,
the chosen atom is isolated onto its own species label (same Z, own
pseudopotential file) behind the scenes; --species stays the real chemical
species everywhere else (the final recommended %%block LDAU.proj, the
REFERENCE_U sanity check, ...).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --species Mn\n"
               "  %(prog)s --species Fe --shell 3d --j 0.5\n"
               "  %(prog)s --species Mn --pseudo-dir ./pseudos\n"
               "  %(prog)s --species Mn --pseudo-dir dojo --save-report\n"
               "  %(prog)s --species Mn --atom-index 1  "
               "# Mn appears twice in structure.fdf\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure.fdf, copied as 'structure.fdf' into every run "
                             "folder (default: structure.fdf).")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                        help="Existing calc.fdf template (default: calc.fdf), with %%include "
                             "structure.fdf and everything else already configured; the DFT+U "
                             "perturbation lines (and single-point SCF enforcement) are added "
                             "via a config_extra.fdf %%include, not embedded inline.")
    parser.add_argument("--species", type=str, required=True,
                        help="Species label to apply the Hubbard U correction to (must be present in "
                             "--structure).")
    parser.add_argument("--atom-index", type=int, default=None, metavar="N",
                        help="1-based index (counting ALL atoms in --structure's coordinate block, in "
                             "file order -- same convention as stb-defect --index) of the specific atom "
                             "to perturb. Required only if --species appears more than once in "
                             "--structure; omit it when --species matches a single atom.")
    parser.add_argument("--shell", type=str, default=None, choices=sorted(SHELL_NAMES),
                        help="Correlated shell (3d/4d/5d/4f/5f). Default: the standard shell for "
                             "--species (transition metals -> (n)d, lanthanides -> 4f, actinides -> 5f).")
    parser.add_argument("--j", type=float, default=0.0,
                         help="Exchange J (eV), fixed across the whole workflow (default: 0.0).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="", metavar="DIR",
                         help="Where to get the pseudopotential files (.psf/.psml) SIESTA needs -- "
                              f"a bundled bank ({', '.join(BANKS)}) or a folder path. Only the "
                              "species actually present in --structure are copied into 'reference/' "
                              "now, and into every scf_alpha_*/frozen_alpha_* folder later by "
                              "stb-hubbarduAlphas. If omitted, you must copy them into every run "
                              "folder yourself.")
    parser.add_argument("-o", "--output-dir", type=str, default="hubbardu_runs",
                        help="Output directory for the run folders (default: hubbardu_runs).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbardu {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    # os.path.isdir/exists never expand '~' themselves -- that's a shell-only
    # behavior, and argv passed through subprocess.run() (as the interactive
    # stb-suite menu does) never goes through a shell at all, so a path like
    # '~/pseudos' would otherwise be checked completely literally and always
    # fail.
    args.structure = os.path.expanduser(args.structure)
    args.calc = os.path.expanduser(args.calc)
    args.output_dir = os.path.expanduser(args.output_dir)

    if not os.path.exists(args.structure):
        sys.exit(color_text(f"Error: File '{args.structure}' not found.", 'red'))
    if not os.path.exists(args.calc):
        sys.exit(color_text(f"Error: File '{args.calc}' not found.", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"Error: {e}", 'red'))

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-HUBBARDU STAGE 1 REPORT (HUBBARD U REFERENCE PREP) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file         : {args.structure}", f_out)
    print_dual(f"Calc.fdf template      : {args.calc}", f_out)
    print_dual(f"Output directory       : {args.output_dir}", f_out)
    print_dual(f"Species                : {args.species}", f_out)
    print_dual(f"Atom index             : "
               f"{args.atom_index if args.atom_index is not None else 'auto (single atom)'}", f_out)
    print_dual(f"Shell                  : {args.shell or 'auto-detect'}", f_out)
    print_dual(f"Exchange J (eV)        : {args.j:.3f}", f_out)
    print_dual(f"Pseudopotentials       : {args.pseudo_dir or 'not given'}", f_out)
    print_dual(f"Report                 : {report_path if report_path else '(not saved)'}", f_out)

    print_section("[1] ATOM SELECTION", f_out)
    try:
        structure = structure_io.read_fdf(args.structure)
    except (FileNotFoundError, ValueError) as e:
        fail(str(e))

    if args.species not in structure.species:
        fail(f"Species '{args.species}' not found in '{args.structure}' "
             f"(available: {', '.join(structure.species)}).")

    n_atoms = len(structure.atoms)
    species_atom_indices = [i for i, (symbol, _) in enumerate(structure.atoms) if symbol == args.species]
    alias_label = None
    print_dual(f"Total atoms in structure : {n_atoms}", f_out)
    print_dual(f"Atom(s) of '{args.species}' : {len(species_atom_indices)} found "
               f"(index(es) {', '.join(str(i + 1) for i in species_atom_indices)})", f_out)

    if len(species_atom_indices) > 1:
        if args.atom_index is None:
            print_dual(color_text(
                f"[FAIL] --species '{args.species}' appears {len(species_atom_indices)} times in "
                f"'{args.structure}' -- the linear-response perturbation must be isolated to a single "
                "atom (perturbing every atom of the species at once measures a collective response, "
                "not the isolated single-site response the method requires). Pass --atom-index to "
                "pick one:", 'red'), f_out)
            try:
                pmg_structure = structure_io.to_pymatgen(structure)
                sites, space_group = find_inequivalent_sites(pmg_structure, symprec=0.01,
                                                               filter_species=args.species)
                print_dual(f"Space group: {space_group}", f_out)
                print_table(["--atom-index", "Wyckoff", "Note"], [
                    ([str(rep + 1), wyckoff,
                      (f"{multiplicity} symmetry-equivalent atom(s) -- any one of them gives the "
                       "same result" if multiplicity > 1 else "symmetrically distinct site")], None)
                    for rep, wyckoff, multiplicity in sites
                ], f_out)
            except ValueError:
                print_table(["--atom-index", "Position"], [
                    ([str(i + 1), str(structure.atoms[i][1])], None) for i in species_atom_indices
                ], f_out)
            if f_out:
                f_out.close()
            sys.exit(1)

        if not (1 <= args.atom_index <= n_atoms):
            fail(f"--atom-index {args.atom_index} is out of range (1 to {n_atoms}).")
        if structure.atoms[args.atom_index - 1][0] != args.species:
            fail(f"atom {args.atom_index} in '{args.structure}' is species "
                 f"'{structure.atoms[args.atom_index - 1][0]}', not '{args.species}'.")

        alias_label = f"{args.species}_pert"
        if alias_label in structure.species:
            fail(f"alias species '{alias_label}' already exists in '{args.structure}' -- rename "
                 "it before running stb-hubbardu on a multi-atom species.")
    elif args.atom_index is not None:
        if not (1 <= args.atom_index <= n_atoms) or structure.atoms[args.atom_index - 1][0] != args.species:
            fail(f"--atom-index {args.atom_index} does not refer to species '{args.species}' "
                 f"in '{args.structure}'.")

    perturbed_species = alias_label or args.species

    if alias_label:
        print_dual(f"Perturbed atom         : atom #{args.atom_index} of '{args.species}', isolated "
                   f"via alias species '{alias_label}' (same Z, own pseudopotential) -- the other "
                   f"atom(s) of '{args.species}' are left unperturbed.", f_out)
        print_dual(color_text(
            f"[NOTE] If your calc.fdf template configures a custom %block PAO.Basis/PAO.Basis.Sizes "
            f"for '{args.species}', it is NOT automatically applied to '{alias_label}' -- SIESTA will "
            "auto-generate a basis for it instead. Check convergence/consistency if you rely on "
            "non-default basis parameters.", 'yellow'), f_out)
    else:
        print_dual(f"Perturbing the only atom of '{args.species}' "
                   f"(index {species_atom_indices[0] + 1}) -- no alias needed.", f_out)

    print_section("[2] DFT+U PERTURBATION SETUP", f_out)
    shell_name = args.shell or DEFAULT_SHELL.get(args.species)
    if shell_name is None:
        fail(f"No default correlated shell known for '{args.species}' -- "
             f"pass --shell explicitly ({', '.join(sorted(SHELL_NAMES))}).")
    n, l = SHELL_NAMES[shell_name]

    with open(args.calc, 'r') as f:
        calc_template = f.read()
    check_no_existing_ldau(calc_template, args.calc)
    with open(args.structure, 'r') as f:
        structure_text = f.read()
    label = get_system_label(calc_template, structure_text)

    print_dual(f"SystemLabel            : {label}", f_out)
    print_dual(f"Species / shell        : {args.species} ({shell_name}: n={n}, l={l})", f_out)
    print_dual(f"Perturbation           : alpha=0.0 (reference run)", f_out)
    print_dual(f"Forced to single-point SCF via '%include {EXTRA_FDF_FILE}'.", f_out)

    ref_extra = ldau_perturbation_block(perturbed_species, n, l, 0.0, args.j, potential_shift=True)
    preview_extra = ("# Auto-generated by stb-hubbardu -- DFT+U perturbation "
                      "(Cococcioni-de Gironcoli) + single-point SCF enforcement.\n"
                      + ref_extra + "\n" + SINGLE_POINT_BLOCK)
    print_dual(f"\n{EXTRA_FDF_FILE} (written into every generated folder):", f_out)
    for line in preview_extra.rstrip("\n").split("\n"):
        print_dual(f"  {line}", f_out)

    print_section("[3] REFERENCE FOLDER & TEMPLATE SNAPSHOT", f_out)
    # Snapshot the exact template stage 2 must regenerate perturbed folders
    # from -- stage 2 is a separate invocation and must not risk using a
    # different (possibly since-edited) calc.fdf/structure.fdf than the one
    # 'reference/' was actually converged with. When an alias is in play,
    # this snapshot IS the aliased structure, so stage 2 propagates it to
    # every scf_alpha_*/frozen_alpha_* folder for free.
    template_structure_path = os.path.join(args.output_dir, TEMPLATE_STRUCTURE_NAME)
    if alias_label:
        structure_io.alias_single_atom_species(
            args.structure, template_structure_path, args.species, args.atom_index, alias_label)
    else:
        shutil.copy(args.structure, template_structure_path)
    with open(os.path.join(args.output_dir, TEMPLATE_CALC_NAME), 'w') as f:
        f.write(calc_template)
    print_dual(color_text(f"[OK] Template snapshot saved: {TEMPLATE_STRUCTURE_NAME}, "
                          f"{TEMPLATE_CALC_NAME}", 'green'), f_out)

    manifest = {
        "species": args.species, "n": n, "l": l, "j": args.j, "label": label,
        "perturbed_species": perturbed_species, "atom_index": args.atom_index,
        "runs": {"reference": {"kind": "reference", "alpha": 0.0}},
    }
    manifest_path = os.path.join(args.output_dir, "run_manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print_dual(color_text(f"[OK] Wrote {manifest_path}", 'green'), f_out)

    ref_folder, _ = write_run_folder(args.output_dir, "reference", calc_template,
                                      template_structure_path, ref_extra)
    print_dual(color_text(
        f"[OK] {ref_folder}/ (structure.fdf, calc.fdf, {EXTRA_FDF_FILE}) -- alpha=0, reference",
        'green'), f_out)

    print_section("[4] PSEUDOPOTENTIALS", f_out)
    if args.pseudo_dir:
        pseudo_snapshot = os.path.join(args.output_dir, PSEUDO_SNAPSHOT_DIRNAME)
        os.makedirs(pseudo_snapshot, exist_ok=True)
        copy_pseudopotentials(args.pseudo_dir, pseudo_snapshot, species=structure.species)
        copied = copy_pseudopotentials(args.pseudo_dir, ref_folder, species=structure.species)
        if alias_label:
            alias_pseudopotential_files(pseudo_snapshot, args.species, alias_label)
            copied = copied + alias_pseudopotential_files(ref_folder, args.species, alias_label)
        if copied:
            print_dual(color_text(
                f"[OK] Copied {len(copied)} pseudopotential file(s) from '{args.pseudo_dir}' into "
                f"'{ref_folder}' (and saved for stb-hubbarduAlphas to reuse)", 'green'), f_out)
        else:
            print_dual(color_text(
                f"[WARNING] --pseudo-dir '{args.pseudo_dir}' has no recognized pseudopotential "
                f"files ({'/'.join(PSEUDO_EXTENSIONS)}) -- did you point it at the right folder? "
                "Nothing was copied.", 'yellow'), f_out)
    elif alias_label:
        print_dual(color_text(
            f"[NOTE] No --pseudo-dir given -- copy your pseudopotential files (.psf/.psml) into "
            "'reference/' yourself before running SIESTA there, including a copy renamed to "
            f"'{alias_label}.<ext>' for the perturbed atom (and again into every scf_alpha_*/"
            "frozen_alpha_* folder stb-hubbarduAlphas generates).", 'yellow'), f_out)
    else:
        print_dual(color_text(
            "[NOTE] No --pseudo-dir given -- copy your pseudopotential files (.psf/.psml) into "
            "'reference/' yourself before running SIESTA there.", 'yellow'), f_out)

    print_section("[5] SUMMARY & NEXT STEPS", f_out)
    print_dual(color_text(f"Success: Generated 'reference/' in '{args.output_dir}'", 'green'), f_out)
    if report_path:
        print_dual(f"Report                 : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in '{args.output_dir}/reference/'.", f_out)
    print_dual(f"  2. Once it converges, run: stb-hubbarduAlphas --dir {args.output_dir}", f_out)
    print_dual(f"     (it will read '{label}.DM' from reference/ and copy it in for you)", f_out)

    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
