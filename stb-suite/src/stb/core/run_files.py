"""What the files of a SIESTA run folder are, and which files its .fdf input reads.

Shared by stb-clean (3.2: what stays and what goes) and stb-archive (3.7: what goes into the archive) -- moved here from
clean.py when stb-archive became the second consumer (the suite's extract-on-second-use rule). The kind of a file
comes from SIESTA's own naming conventions; read_fdf_settings follows %include, '%block X < file' and 'Label < file'
the way SIESTA does (first occurrence of a label wins).
"""

import os
import re


CAT_RESTART = "restart (SIESTA reads it back)"
CAT_LOG = "SIESTA log / status"
CAT_JOB = "job-scheduler output"
CAT_GEOM = "geometry / forces / dynamics"
CAT_ELEC = "electronic structure"
CAT_GRID = "real-space grids"
CAT_BASIS = "basis by-products"
CAT_POST = "post-processing / plots"
CAT_SCRIPT = "job script (not in --keep)"
CAT_OTHER = "other"
CATEGORY_ORDER = [CAT_RESTART, CAT_LOG, CAT_JOB, CAT_GEOM, CAT_ELEC, CAT_GRID, CAT_BASIS, CAT_POST, CAT_SCRIPT,
                  CAT_OTHER]

_RESTART_EXT = {".DM", ".XV", ".CG", ".LWF", ".TSDE", ".TSHS", ".ZM", ".STRUCT_NEXT_ITER", ".MDX"}
_LOG_NAMES = {"MESSAGES", "TIMES", "CLOCK", "0_NORMAL_EXIT", "PARALLEL_DIST", "NON_TRIMMED_KP_LIST",
              "OUTVARS.yml", "FORCE_STRESS", "BASIS_ENTHALPY", "BASIS_HARRIS_ENTHALPY", "INPUT_TMP"}
_LOG_EXT = {".out", ".log", ".times", ".alloc", ".bib", ".BASIS_ENTHALPY", ".BASIS_HARRIS_ENTHALPY", ".yml"}
_GEOM_EXT = {".STRUCT_OUT", ".STRUCT_IN", ".ANI", ".MD", ".MDE", ".xyz", ".BONDS", ".BONDS_FINAL", ".FA", ".FAC",
             ".FC", ".ZM_OUT", ".XV_OUT"}
_ELEC_EXT = {".EIG", ".KP", ".bands", ".DOS", ".PDOS", ".ORB_INDX", ".HSX", ".HS", ".WFSX", ".WFS", ".EPSIMG",
             ".BANDS", ".fullBZ", ".PDOS.xml"}
_GRID_EXT = {".RHO", ".VT", ".VH", ".DRHO", ".RHOXC", ".BADER", ".TOCH", ".LDOS", ".VNA", ".RHOINIT", ".IOCH",
             ".RHOINIT", ".nc", ".grid", ".STM", ".CH"}
_POST_EXT = {".dat", ".gplot", ".gnu", ".pdf", ".png", ".svg", ".eps", ".cube", ".xsf", ".ham", ".basis", ".csv",
             ".json", ".txt", ".agr", ".jpg"}
_JOB_RE = re.compile(r"^(slurm-\d+\.out|.*\.[eo]\d+|job\.(o|e|err|out|log)|.*\.err|nohup\.out)$")
_SCRIPT_MARK = re.compile(r"^#(SBATCH|PBS|BSUB|\$ )", re.M)


def classify(name):
    """The kind of file `name` (a basename) is, from SIESTA's own naming conventions."""
    ext = os.path.splitext(name)[1]
    if name.endswith((".ion", ".ion.xml", ".ion.nc")) or name.startswith("basis_set"):
        return CAT_BASIS
    if _JOB_RE.match(name):
        return CAT_JOB
    if ext in _RESTART_EXT or name.endswith(".DM.nc"):
        return CAT_RESTART
    if name in _LOG_NAMES or ext in _LOG_EXT or name.startswith("fdf.") or name.startswith("fdf-"):
        return CAT_LOG
    if ext in _GEOM_EXT:
        return CAT_GEOM
    if ext in _ELEC_EXT or name.endswith(".PDOS.xml"):
        return CAT_ELEC
    if ext in _GRID_EXT:
        return CAT_GRID
    if ext in _POST_EXT:
        return CAT_POST
    return CAT_OTHER


def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} GB"


def file_size(path):
    try:
        return os.lstat(path).st_size
    except OSError:
        return 0


def looks_like_job_script(path):
    """True for a small text file that starts with '#!' or carries scheduler directives / an mpirun line."""
    try:
        if os.path.getsize(path) > 1_000_000:
            return False
        with open(path, "r", errors="strict") as f:
            head = f.read(4000)
    except (OSError, UnicodeDecodeError):
        return False
    return head.startswith("#!") or bool(_SCRIPT_MARK.search(head))


# ---------------------------------------------------------------- the fdf input ---------------------------------
def fdf_key(label):
    return re.sub(r"[._-]", "", label).lower()


def fdf_true(value):
    return value is None or value.strip().lower().strip(".") in ("t", "true", "yes", "y")


def read_fdf_settings(fdf_paths):
    """Reads the kept .fdf files (in the order given, first occurrence wins, like SIESTA) for what matters to a
    restart: the files they read, SystemLabel, the UseSave* flags and the species labels."""
    refs, labels, species = {}, {}, []
    in_species = False
    for path in fdf_paths:
        try:
            lines = open(path, errors="replace").read().splitlines()
        except OSError:
            continue
        base = os.path.dirname(path)
        for raw in lines:
            line = re.split(r"[#!;]", raw, maxsplit=1)[0].strip()
            if not line:
                continue
            low = line.lower()
            m = re.match(r"%include\s+(\S+)", line, re.I)
            if m:
                refs.setdefault(os.path.normpath(os.path.join(base, m.group(1))), path)
                continue
            m = re.search(r"<\s*(\S+)\s*$", line)
            if m:
                refs.setdefault(os.path.normpath(os.path.join(base, m.group(1))), path)
            if low.startswith("%block"):
                in_species = fdf_key(line.split()[1]) == "chemicalspecieslabel" if len(line.split()) > 1 else False
                continue
            if low.startswith("%endblock"):
                in_species = False
                continue
            if in_species:
                parts = line.split()
                if len(parts) >= 3:
                    species.append(parts[2])
                continue
            parts = line.split(None, 1)
            labels.setdefault(fdf_key(parts[0]), parts[1] if len(parts) > 1 else None)
    return refs, labels, list(dict.fromkeys(species))



def collect_files(path, recursive):
    """Every regular file (or symlink to one) under `path`; hidden directories are never entered."""
    found = []
    if recursive:
        for root, dirs, files in os.walk(path):
            dirs[:] = sorted(d for d in dirs if not d.startswith('.'))
            found += [os.path.join(root, f) for f in files]
    else:
        found = [os.path.join(path, f) for f in os.listdir(path)]
    return sorted(f for f in found if os.path.isfile(f))
