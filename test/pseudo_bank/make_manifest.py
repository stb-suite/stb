#!/usr/bin/env python
"""Writes manifest.json in every folder of the bundled PseudoDojo bank (maintainer tool).

Run after the bank changes:
    python test/pseudo_bank/make_manifest.py --archives DIR

DIR holds the official archives (http://www.pseudo-dojo.org/pseudos/<table>_psml.tgz). Each
folder's files are compared byte by byte with its archive, and the archive's SHA-256 is recorded.
Without --archives the existing 'source' entries are kept and only the per-file data is rebuilt.
"""
import argparse
import hashlib
import json
import os
import sys
import tarfile

from stb.core import pseudopotentials as pp

URL = "http://www.pseudo-dojo.org/pseudos/{}_psml.tgz"
# folder -> (archive stem, date the files were taken from it)
SOURCES = {
    "sr_pbe_0.4_standard": ("nc-sr-04_pbe_standard", "2026-11-20"),
    "sr_pbe_0.4_stringent": ("nc-sr-04_pbe_stringent", "2026-11-20"),
    "sr_pbesol_0.4_standard": ("nc-sr-04_pbesol_standard", "2026-11-17"),
    "sr_pbesol_0.4_stringent": ("nc-sr-04_pbesol_stringent", "2026-11-17"),
    "sr_pw_0.4_standard": ("nc-sr-04_pw_standard", "2026-11-17"),
    "sr_pw_0.4_stringent": ("nc-sr-04_pw_stringent", "2026-11-17"),
    "fr_pbe_0.4_standard": ("nc-fr-04_pbe_standard", "2026-11-17"),
    "fr_pbe_0.4_stringent": ("nc-fr-04_pbe_stringent", "2026-11-17"),
    "fr_pbesol_0.4_standard": ("nc-fr-04_pbesol_standard", "2026-11-17"),
    "fr_pbesol_0.4_stringent": ("nc-fr-04_pbesol_stringent", "2026-11-17"),
    "sr3plus_pbe_0.4_standard": ("nc-sr-04-3plus_pbe_standard", "2026-11-17"),
}


def md5(data):
    return hashlib.md5(data).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archives", help="folder with the official <table>_psml.tgz archives")
    ap.add_argument("--verified", default=None, help="date of the byte-by-byte check (YYYY-MM-DD)")
    args = ap.parse_args()
    root = os.path.dirname(pp.bank_path("dojo"))
    tables = {folder: table for folder, table in pp.DOJO_TABLES.values()}
    failed = False
    for folder in sorted(os.listdir(root)):
        path = os.path.join(root, folder)
        if not os.path.isdir(path):
            continue
        manifest_path = os.path.join(path, pp.BANK_MANIFEST)
        old = json.load(open(manifest_path)) if os.path.isfile(manifest_path) else {}
        stem, taken = SOURCES[folder]
        source = old.get("source", {"archive": stem + "_psml.tgz", "url": URL.format(stem), "taken": taken})
        files = {}
        names = sorted(n for n in os.listdir(path) if n.endswith(".psml"))
        for name in names:
            data = open(os.path.join(path, name), "rb").read()
            h = pp.read_pseudo_header(os.path.join(path, name))
            files[name] = {"md5": md5(data), "z_valence": h["z_valence"], "relativity": h["relativity"],
                           "lj": h["lj"], "scalar_part": h["sr_part"], "ekb_zero": h["n_ekb_zero"]}
        if args.archives:
            tgz = os.path.join(args.archives, source["archive"])
            raw = open(tgz, "rb").read()
            with tarfile.open(tgz) as tar:
                members = {os.path.basename(m.name): m for m in tar.getmembers() if m.name.endswith(".psml")}
                for name in names:
                    if name not in members or md5(tar.extractfile(members[name]).read()) != files[name]["md5"]:
                        print(f"[ERROR] {folder}/{name} differs from {source['archive']}")
                        failed = True
            source["archive_sha256"] = hashlib.sha256(raw).hexdigest()
            source["archive_files"] = len(members)
            if args.verified:
                source["verified"] = args.verified
        manifest = {"schema": "stb_pseudo_bank_v1", "folder": folder, "table": tables[folder],
                    "source": source, "n_files": len(files), "files": files}
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=1, sort_keys=False)
            f.write("\n")
        print(f"{folder}: {len(files)} files, table {tables[folder]}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
