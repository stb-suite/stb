"""Shared atom-identity manifest, consumed by stb-neb.

Meant to be written once by any tool that builds two structures guaranteed
to share atom-index correspondence (one identical-shape copy inside EACH
of the two endpoint folders -- stb-neb's --initial/--final are each given
only ONE folder individually, so a manifest living only at a shared root
would not be discoverable from either folder alone), read and
cross-validated by stb-neb before it trusts that distance-based
(--autosort-tol) atom matching is unnecessary. No in-repo tool currently
produces one (the symmetry-site-enumeration tool that used to, stb-nebSites,
was removed from the NEB workflow), but the format stays as
forward-looking infrastructure for a future producer or a hand-built pair.

Mirrors core/dftu_data.py's REQUIRED_MANIFEST_FIELDS + load_manifest()
"stage 1 writes JSON, stage 2 (a separate, possibly much-later CLI
invocation) loads + validates before trusting anything" pattern -- reused
here since it already solves the identical problem shape for
stb-hubbardu/stb-hubbardUAlphas' own run_manifest.json.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone

MANIFEST_FILENAME = "neb_manifest.json"

# Fields every neb_manifest.json must have -- checked by load_manifest() so
# a malformed/hand-edited/incomplete manifest fails with a clean message
# instead of a raw KeyError somewhere downstream.
REQUIRED_MANIFEST_FIELDS = (
    "schema_version", "pair_id", "site_label", "species_sequence",
    "n_atoms", "n_substrate", "n_adsorbate", "adsorbate", "created",
)

SCHEMA_VERSION = 1


def write_manifest(site_dir, *, pair_id, site_label, species_sequence,
                    n_substrate, n_adsorbate, adsorbate):
    """Writes site_dir/neb_manifest.json.

    `species_sequence` must be the ACTUAL on-disk per-index atom order (the
    caller should re-read it back from the just-written structure.fdf, not
    hand-derive it from whatever pymatgen Structure was written FROM --
    write_fdf() regroups atoms by species before writing, so the physical
    file order is not simply "substrate atoms then adsorbate atoms" in
    general). Returns the manifest path.
    """
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "pair_id": pair_id,
        "site_label": site_label,
        "species_sequence": list(species_sequence),
        "n_atoms": n_substrate + n_adsorbate,
        "n_substrate": n_substrate,
        "n_adsorbate": n_adsorbate,
        "adsorbate": adsorbate,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = os.path.join(site_dir, MANIFEST_FILENAME)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest_path


def load_manifest(manifest_path):
    """Loads and validates a neb_manifest.json, raising ValueError with a
    clean, user-facing message (never a raw KeyError/JSONDecodeError
    surfacing later at some unrelated line) if it's malformed or missing a
    required top-level field -- e.g. a manifest left over from an
    older/incompatible version of this workflow, or hand-edited.
    """
    try:
        with open(manifest_path, "r") as f:
            manifest = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"'{manifest_path}' is not valid JSON ({e}).")
    missing = [field for field in REQUIRED_MANIFEST_FIELDS if field not in manifest]
    if missing:
        raise ValueError(
            f"'{manifest_path}' is missing required field(s) {missing} -- this manifest is "
            "either from an older/incompatible schema version or a hand-edited file."
        )
    return manifest


def validate_manifest_pair(manifest_initial, manifest_final):
    """Cross-checks two already-loaded manifests against EACH OTHER.

    The hard guarantee this proves is `species_sequence` being LITERALLY
    identical (list equality, order matters) between the two -- that's the
    guarantee any producer of a matched manifest pair (e.g. two structures
    built from the same base structure, only the endpoint geometry
    differing) needs to provide, and it alone is sufficient proof of a
    correct index correspondence. Raises ValueError if the sequences differ
    in length or content (reports up to the first 10 mismatched indices).

    Returns True/False for whether `pair_id` ALSO matches -- a softer,
    purely-informational "these two folders came from the same producing
    run" signal. A pair_id mismatch does NOT raise: species_sequence
    equality alone is sufficient proof of order, regardless of which run
    produced either folder.
    """
    seq_i = manifest_initial["species_sequence"]
    seq_f = manifest_final["species_sequence"]
    if len(seq_i) != len(seq_f):
        raise ValueError(
            f"--initial's manifest has {len(seq_i)} atom(s), --final's has {len(seq_f)} -- "
            "these do not look like a matching manifest pair."
        )
    mismatches = [i for i, (a, b) in enumerate(zip(seq_i, seq_f)) if a != b]
    if mismatches:
        shown = mismatches[:10]
        raise ValueError(
            f"--initial/--final manifests disagree on species_sequence at index(es) "
            f"{shown}{'...' if len(mismatches) > 10 else ''} -- these do not look like a "
            "matching manifest pair (or one of the two structure.fdf files was "
            "regenerated/hand-edited since the manifest was written)."
        )
    return manifest_initial["pair_id"] == manifest_final["pair_id"]
