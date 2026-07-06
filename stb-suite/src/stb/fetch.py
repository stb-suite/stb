#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.13.0"

import sys
import argparse
import requests
from pymatgen.core import Composition, Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, get_input, show_intro
from stb.core.symmetry import UNITCELL_MODES, reduce_to_unitcell

COD_BASE_URL = "https://www.crystallography.net/cod"

# A hand-picked subset of pymatgen's ~30 built-in OptimadeRester aliases,
# relevant to a DFT/SIESTA audience (mp/cod deliberately excluded -- they're
# already better served by this tool's own dedicated --source materials
# -project/cod integrations, e.g. energy_above_hull-based --most-stable
# ranking that generic OPTIMADE fields don't cleanly give us). OPTIMADE
# endpoints are third-party services that move/break independently of
# pymatgen's release cycle -- "confirmed working"/"currently unreliable"
# notes below reflect a live check done when this list was written, not a
# permanent guarantee.
CURATED_OPTIMADE_PROVIDERS = [
    ("twodmatpedia", "2D Materials Encyclopedia (2D materials) -- confirmed working"),
    ("jarvis", "NIST JARVIS-DFT (2D + 3D materials)"),
    ("aflow", "AFLOW"),
    ("oqmd", "Open Quantum Materials Database"),
    ("alexandria.pbe", "Alexandria (PBE)"),
    ("odbx", "Open Database of Xtals"),
    ("mcloud.mc2d", "Materials Cloud 2D Structures Database"),
    ("mcloud.2dtopo", "Materials Cloud 2D Topological Database"),
    ("cmr", "DTU Computational Materials Repository, incl. C2DB -- currently unreliable, endpoint may be stale"),
]


def to_hill_formula(formula):
    """COD's and OPTIMADE's formula search want element symbols in Hill
    notation separated by spaces (e.g. 'Fe3 O4'). Accepts a formula typed
    any way pymatgen's Composition can parse ('Fe3O4', 'Fe3 O4', ...) and
    reformats it to match.
    """
    return Composition(formula).hill_formula


def resolve_disorder_or_fail(structure, source, entry_id):
    """Real database entries (COD especially) are often disordered even for
    common materials -- e.g. magnetite's octahedral Fe site is frequently
    given as a Fe2+/Fe3+ mixed-occupancy split rather than plain Fe, since
    that's the crystallographically precise description. That kind of
    same-element-different-oxidation-state disorder is safe to collapse to a
    single element at full occupancy for a plain (non charge-ordered) SIESTA
    structure. Genuine chemical disorder (different elements sharing a site)
    can't be resolved automatically -- raises ValueError in that case.

    Returns (structure, was_collapsed).
    """
    if structure.is_ordered:
        return structure, False

    species = []
    for site in structure:
        symbols = {sp.symbol for sp in site.species}
        if len(symbols) != 1:
            raise ValueError(
                f"{source} entry {entry_id} has genuine chemical disorder (a site mixes "
                f"{', '.join(sorted(symbols))}) -- this can't be written as a definite-atom "
                ".fdf. Try a different/ordered entry."
            )
        species.append(symbols.pop())

    ordered = Structure(structure.lattice, species, structure.frac_coords)
    return ordered, True


def fetch_cod_candidates(formula):
    """Returns the full list of rows from COD's formula search (no capping --
    that's a display concern for the caller), where each row is a dict with
    the fields COD returns (file, formula, sg, sgNumber, a/b/c, year...).
    """
    response = requests.get(
        f"{COD_BASE_URL}/result",
        params={"formula": to_hill_formula(formula), "format": "json"},
        timeout=30,
    )
    response.raise_for_status()
    return response.json()


def fetch_cod_structure(cod_id):
    """Downloads and parses a single COD entry's CIF into a pymatgen Structure."""
    response = requests.get(f"{COD_BASE_URL}/{cod_id}.cif", timeout=30)
    response.raise_for_status()
    structures = CifParser.from_str(response.text).parse_structures(primitive=False, on_error="warn")
    if not structures:
        raise ValueError(f"COD entry {cod_id} did not contain a parseable structure.")
    return structures[0]


def fetch_mp_candidates(mpr, formula):
    """Returns the full list of rows from Materials Project's summary search
    (no capping -- that's a display concern for the caller)."""
    fields = ["material_id", "formula_pretty", "symmetry", "energy_above_hull", "nsites"]
    return mpr.get_summary({"formula": formula}, fields=fields)


def _single_provider_structures(results):
    """OptimadeRester results are keyed by provider name, but since we always
    query exactly one provider there's at most one key -- return its
    {id: Structure} dict (or {} if the provider returned nothing/failed).
    """
    return next(iter(results.values()), {})


def fetch_optimade_candidates(provider, formula):
    """Returns a list of {"id", "formula", "natoms", "_structure"} dicts from
    an OPTIMADE provider's formula search. Unlike COD/MP, OPTIMADE already
    returns parsed Structure objects at search time; stashing them under
    "_structure" lets picking a candidate skip a second network round-trip.

    Searches by element set + element count (elements/nelements) rather than
    the OPTIMADE chemical_formula_hill filter: confirmed live that at least
    one real provider (twodmatpedia) errors on chemical_formula_hill queries
    (a server-side quirk, not universally implemented the same way), while
    elements/nelements is a core OPTIMADE filter field every provider must
    support. Trade-off: this matches every stoichiometry of the same
    elements (e.g. searching 'MoS2' also returns MoS, Mo3S4, ...), not just
    the exact formula -- the candidate table shows each result's actual
    formula so the user can still pick the right one.
    """
    from pymatgen.ext.optimade import OptimadeRester

    elements = [str(el) for el in Composition(formula).elements]
    with OptimadeRester(provider) as opt:
        results = opt.get_structures(elements=elements, nelements=len(elements))
    structures = _single_provider_structures(results)
    return [
        {"id": entry_id, "formula": s.composition.reduced_formula, "natoms": len(s), "_structure": s}
        for entry_id, s in structures.items()
    ]


def fetch_optimade_structure(provider, optimade_id):
    """Fetches a single OPTIMADE entry by its provider-specific id."""
    from pymatgen.ext.optimade import OptimadeRester

    with OptimadeRester(provider) as opt:
        results = opt.get_structures_with_filter(f'id="{optimade_id}"')
    structures = _single_provider_structures(results)
    if not structures:
        raise ValueError(f"OPTIMADE provider '{provider}' has no entry with id '{optimade_id}'.")
    return next(iter(structures.values()))


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Fetches a structure from an online database and writes it as .fdf.", 'bold')}
Look it up by exact id, or search by formula and pick a candidate.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --source cod --cod-id 1010369\n"
               "  %(prog)s --source cod --formula 'Fe3O4'\n"
               "  %(prog)s --source materials-project --material-id mp-19306\n"
               "  %(prog)s --source materials-project --formula Fe3O4 --most-stable\n"
               "  %(prog)s --source optimade --provider twodmatpedia --formula MoS2\n"
               "  %(prog)s --source optimade --provider twodmatpedia --optimade-id 2dm-2127\n"
               "  %(prog)s --source cod --cod-id 1010369 --unitcell primitive\n"
               "  %(prog)s --list-providers\n"
    )

    parser.add_argument("--source", choices=["materials-project", "cod", "optimade"], default=None,
                        help="Database to query. Required unless --list-providers is given.")

    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--material-id", type=str, default=None,
                           help="Exact Materials Project id (e.g. 'mp-19306'). materials-project only.")
    selection.add_argument("--cod-id", type=int, default=None,
                           help="Exact Crystallography Open Database numeric id. cod only.")
    selection.add_argument("--optimade-id", type=str, default=None,
                           help="Exact id within the chosen --provider. optimade only.")
    selection.add_argument("--formula", type=str, default=None,
                           help="Search by chemical formula (e.g. 'Fe3O4'). Works with any source; "
                                "prints a candidate table if more than one match is found.")

    parser.add_argument("--most-stable", action="store_true",
                        help="With --formula on materials-project, auto-pick the candidate with the "
                             "lowest energy above hull instead of requiring an explicit --material-id. "
                             "Not available for --source cod/optimade (no DFT energy to rank by).")
    parser.add_argument("--api-key", type=str, default=None,
                        help="Materials Project API key. Default: the PMG_MAPI_KEY environment "
                             "variable / pymatgen's own settings file.")
    parser.add_argument("--provider", type=str, default=None,
                        help="OPTIMADE provider: an alias (see --list-providers for a curated "
                             "subset; any alias from pymatgen's own OptimadeRester.aliases also "
                             "works) or a raw OPTIMADE base URL. Required with --source optimade.")
    parser.add_argument("--list-providers", action="store_true",
                        help="Print a curated list of OPTIMADE providers and exit. Does not "
                             "require --source or a selection flag.")
    parser.add_argument("--limit", type=int, default=20,
                        help="Cap on how many --formula candidates to display (default: 20).")

    parser.add_argument("--unitcell", choices=list(UNITCELL_MODES), default=None,
                        help="Also reduce the fetched structure to its primitive/conventional/"
                             "refined unit cell (same logic as stb-unitcell) before writing. "
                             "Default: keep the structure exactly as fetched.")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision, used for the fetched structure's space group "
                             "summary and for --unitcell (default: 1e-3, matches stb-unitcell/"
                             "stb-symmetry).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry angle tolerance in degrees, used the same way as "
                             "--symprec (default: 5.0, pymatgen's own default).")

    parser.add_argument("-o", "--output", type=str, default="fetched.fdf",
                        help="Output .fdf file name (default: fetched.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-fetch {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.list_providers:
        print(color_text("Curated OPTIMADE providers:", 'bold'))
        for alias, desc in CURATED_OPTIMADE_PROVIDERS:
            print(f"  {color_text(alias, 'cyan'):<20} {desc}")
        print("\nAny other alias from pymatgen's OptimadeRester.aliases, or a raw OPTIMADE "
              "base URL, also works with --provider.")
        sys.exit(0)

    if args.source is None:
        parser.error("--source is required (unless --list-providers is given).")
    if not any([args.material_id, args.cod_id, args.optimade_id, args.formula]):
        parser.error("one of --material-id, --cod-id, --optimade-id, --formula is required.")

    if args.source == "cod":
        if args.material_id is not None:
            parser.error("--material-id is only valid with --source materials-project.")
        if args.api_key is not None:
            parser.error("--api-key is only valid with --source materials-project.")
        if args.most_stable:
            parser.error("--most-stable is only valid with --source materials-project.")
        if args.provider is not None:
            parser.error("--provider is only valid with --source optimade.")
        if args.optimade_id is not None:
            parser.error("--optimade-id is only valid with --source optimade.")
    elif args.source == "materials-project":
        if args.cod_id is not None:
            parser.error("--cod-id is only valid with --source cod.")
        if args.provider is not None:
            parser.error("--provider is only valid with --source optimade.")
        if args.optimade_id is not None:
            parser.error("--optimade-id is only valid with --source optimade.")
    else:  # optimade
        if args.cod_id is not None:
            parser.error("--cod-id is only valid with --source cod.")
        if args.material_id is not None:
            parser.error("--material-id is only valid with --source materials-project.")
        if args.api_key is not None:
            parser.error("--api-key is only valid with --source materials-project.")
        if args.most_stable:
            parser.error("--most-stable is only valid with --source materials-project.")
        if args.provider is None:
            parser.error("--source optimade requires --provider (see --list-providers).")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Fetch a structure from an online database:", 'bold'))
    print("-" * 60)
    print(f"  {color_text('Source:', 'cyan')} {args.source}")

    material_id = args.material_id
    cod_id = args.cod_id
    optimade_id = args.optimade_id
    pending_structure = None
    mpr = None

    if args.source == "materials-project":
        from pymatgen.ext.matproj import MPRester, MPRestError
        try:
            mpr = MPRester(api_key=args.api_key)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

    if args.formula:
        print(f"  {color_text('Formula search:', 'cyan')} {args.formula}")
        try:
            if args.source == "cod":
                all_rows = fetch_cod_candidates(args.formula)
            elif args.source == "materials-project":
                all_rows = fetch_mp_candidates(mpr, args.formula)
            else:
                all_rows = fetch_optimade_candidates(args.provider, args.formula)
        except requests.exceptions.RequestException as e:
            print(color_text(f"Error: network request failed: {e}", 'red'))
            sys.exit(1)
        except Exception as e:  # pymatgen's MPRestError, or any decode issue
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

        total = len(all_rows)
        if total == 0:
            print(color_text(f"Error: no {args.source} entries found for formula '{args.formula}'.", 'red'))
            sys.exit(1)

        # --most-stable ranks over every match found; the --limit cap below is
        # purely a display/interactive-pick concern, not a search restriction.
        rows = all_rows[:args.limit]

        print(f"\n  {color_text('Found', 'cyan')} {total} candidate(s)"
              f"{' (showing first ' + str(len(rows)) + ')' if total > len(rows) else ''}:")
        if args.source == "cod":
            for i, row in enumerate(rows):
                print(f"    {color_text(f'[{i + 1}]', 'cyan'):<5} COD {row['file']:<10} {row['formula']:<16} "
                      f"sg={row.get('sg', '?'):<14} a={row.get('a', '?')}")
        elif args.source == "materials-project":
            for i, row in enumerate(rows):
                symmetry = row.get("symmetry") or {}
                sg_symbol = symmetry.get("symbol", "?") if isinstance(symmetry, dict) else "?"
                print(f"    {color_text(f'[{i + 1}]', 'cyan'):<5} {row['material_id']:<12} "
                      f"{row.get('formula_pretty', '?'):<16} sg={sg_symbol:<14} "
                      f"E_above_hull={row.get('energy_above_hull', '?')}")
        else:
            for i, row in enumerate(rows):
                print(f"    {color_text(f'[{i + 1}]', 'cyan'):<5} {row['id']:<16} {row['formula']:<16} "
                      f"{row['natoms']} atoms")

        if total > 1:
            if args.source == "materials-project" and args.most_stable:
                best = min(all_rows, key=lambda r: r.get("energy_above_hull", float("inf")))
                material_id = best["material_id"]
                print(f"\n  {color_text('--most-stable selected:', 'yellow')} {material_id}"
                      f" (E_above_hull={best.get('energy_above_hull', '?')}, out of all {total} matches)")
            elif sys.stdin.isatty():
                selected = None
                while selected is None:
                    choice = get_input(f"\nPick a candidate [1-{len(rows)}] (blank to cancel): ").strip()
                    if not choice:
                        print(color_text("Cancelled.", 'red'))
                        sys.exit(1)
                    if choice.isdigit() and 1 <= int(choice) <= len(rows):
                        selected = rows[int(choice) - 1]
                    else:
                        print(color_text(f"Enter a number between 1 and {len(rows)}.", 'red'))
                if args.source == "materials-project":
                    material_id = selected["material_id"]
                elif args.source == "cod":
                    cod_id = int(selected["file"])
                else:
                    optimade_id = selected["id"]
                    pending_structure = selected["_structure"]
            else:
                flag = {"materials-project": "--material-id", "cod": "--cod-id",
                        "optimade": "--optimade-id"}[args.source]
                print(color_text(
                    f"\nError: multiple candidates found; rerun with {flag} to pick one"
                    + (" or add --most-stable." if args.source == "materials-project" else "."), 'red'))
                sys.exit(1)
        else:
            row = rows[0]
            if args.source == "materials-project":
                material_id = row["material_id"]
            elif args.source == "cod":
                cod_id = int(row["file"])
            else:
                optimade_id = row["id"]
                pending_structure = row["_structure"]

    print(f"\n  {color_text('Fetching...', 'yellow')}")
    try:
        if args.source == "materials-project":
            try:
                structure = mpr.get_structure_by_material_id(material_id)
            except (MPRestError, IndexError) as e:
                print(color_text(f"Error: could not fetch '{material_id}': {e}", 'red'))
                sys.exit(1)
        elif args.source == "cod":
            structure = fetch_cod_structure(cod_id)
        else:
            structure = pending_structure if pending_structure is not None \
                else fetch_optimade_structure(args.provider, optimade_id)
    except requests.exceptions.RequestException as e:
        print(color_text(f"Error: network request failed: {e}", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    fetched_id = {"materials-project": material_id, "cod": cod_id, "optimade": optimade_id}[args.source]

    try:
        structure, was_collapsed = resolve_disorder_or_fail(structure, args.source, fetched_id)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)
    if was_collapsed:
        print(color_text(
            "  Note: entry had oxidation-state-only disorder (e.g. Fe2+/Fe3+ on the same "
            "site) -- collapsed to a single element per site.", 'yellow'))

    sga = SpacegroupAnalyzer(structure, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
    print(f"\n  {color_text('Fetched:', 'cyan')} {args.source} {fetched_id}")
    print(f"  {color_text('Formula:', 'cyan')} {structure.composition.reduced_formula}")
    print(f"  {color_text('Atoms:', 'cyan')} {len(structure)}")
    print(f"  {color_text('Space group:', 'cyan')} {sga.get_space_group_symbol()} (No. {sga.get_space_group_number()})")

    if args.unitcell:
        fetched_atoms = len(structure)
        try:
            structure, sga = reduce_to_unitcell(
                structure, args.unitcell, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
        except ValueError as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)
        print(f"\n  {color_text('Unit cell mode:', 'cyan')} {args.unitcell}")
        print(f"  {color_text('Output atoms:', 'cyan')} {len(structure)}")
        if args.unitcell != "refined" and len(structure) == fetched_atoms:
            print(color_text(
                f"  Note: fetched structure is already the {args.unitcell} cell at this "
                "symprec (no reduction).", 'yellow'))

    species_meta = {}
    for symbol in dict.fromkeys(site.specie.symbol for site in structure):
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(structure, species_meta=species_meta)
    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
