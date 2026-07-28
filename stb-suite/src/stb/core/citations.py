"""Shared BibTeX citation entries and writer, used by every stb tool that
writes a references.bib alongside its output (currently stb-inputfile,
stb-kgrid, and stb-kpath). Each entry is a (key, bibtex-entry-string) tuple;
tool-specific entries (DFT-D3, Nose thermostat, pseudopotential banks,
Setyawan-Curtarolo, ...) stay local to their own tool module.
"""

import os
import re

SIESTA = ("Soler2002", """@article{Soler2002,
  author  = {Soler, Jos\\'e M. and Artacho, Emilio and Gale, Julian D. and Garc\\'ia, Alberto and Junquera, Javier and Ordej\\'on, Pablo and S\\'anchez-Portal, Daniel},
  title   = {The {SIESTA} method for ab initio order-{N} materials simulation},
  journal = {Journal of Physics: Condensed Matter},
  year    = {2002},
  volume  = {14},
  number  = {11},
  pages   = {2745--2779},
  doi     = {10.1088/0953-8984/14/11/302}
}""")

SIESTA_RECENT = ("Garcia2020", """@article{Garcia2020,
  author  = {Garc\\'ia, Alberto and Papior, Nick and Akhtar, Arsalan and Nordstr\\"om, Emilio and Ordej\\'on, Pablo and Blum, Volker and Sanchez-Portal, Daniel},
  title   = {Siesta: Recent developments and applications},
  journal = {The Journal of Chemical Physics},
  year    = {2020},
  volume  = {152},
  number  = {20},
  pages   = {204108},
  doi     = {10.1063/5.0005077}
}""")

MONKHORST_PACK = ("MonkhorstPack1976", """@article{MonkhorstPack1976,
  author  = {Monkhorst, Hendrik J. and Pack, James D.},
  title   = {Special points for {Brillouin-zone} integrations},
  journal = {Physical Review B},
  year    = {1976},
  volume  = {13},
  number  = {12},
  pages   = {5188--5192},
  doi     = {10.1103/PhysRevB.13.5188}
}""")

# The MACE architecture itself -- cite regardless of whether the foundation
# checkpoint or a custom fine-tuned model is used (a fine-tuned model is
# still this architecture). Moved here from mlrelax.py once stacking2D.py's
# --ml-relax became a second consumer of the exact same entry.
MACE = ("Batatia2022mace", """@inproceedings{Batatia2022mace,
  author    = {Batatia, Ilyes and Kov\\'acs, D\\'avid P\\'eter and Simm, Gregor N. C. and Ortner, Christoph and Cs\\'anyi, G\\'abor},
  title     = {{MACE}: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields},
  booktitle = {Advances in Neural Information Processing Systems (NeurIPS)},
  year      = {2022},
  eprint    = {2206.07697},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ML}
}""")

# The MACE-MP-0 foundation model/training-data paper -- only cite when the
# foundation checkpoint is actually used (not a --custom-model).
MACE_MP = ("Batatia2023foundation", """@misc{Batatia2023foundation,
  author = {Batatia, Ilyes and Benner, Philipp and Chiang, Yuan and Elena, Alin M. and Kov\\'acs, D\\'avid P\\'eter and Riebesell, Janosh and Advincula, Xavier R. and Asta, Mark and Baldwin, William J. and Bernstein, Noam and Bhowmik, Arghya and Blau, Samuel M. and C\\u{a}rare, Vlad and Darby, James P. and De, Sandip and Della Pia, Flaviano and Deringer, Volker L. and Elijo\\v{s}ius, Rokas and El-Machachi, Zakariya and Fako, Edvin and Ferrari, Andrea C. and Genreith-Schriever, Annalena and George, Janine and Goodall, Rhys E. A. and Grey, Clare P. and Han, Shuang and Handley, Will and Heenen, Hendrik H. and Hermansson, Kersti and Holm, Christian and Jaafar, Jad and Hofmann, Stephan and Jakob, Konstantin S. and Jung, Hyunwook and Kapil, Venkat and Kaplan, Aaron D. and Karimitari, Nima and Kroupa, Namu and Kullgren, Jolla and Kuner, Matthew C. and Kuryla, Domantas and Liepuoniute, Guoda and Margraf, Johannes T. and others},
  title  = {A foundation model for atomistic materials chemistry},
  year   = {2023},
  eprint = {2401.00096},
  archivePrefix = {arXiv},
  primaryClass  = {physics.chem-ph}
}""")

# icet, the cluster-expansion library stb-sqs uses for its SQS search
# (Monte Carlo or enumeration).
ICET = ("Angqvist2019icet", """@article{Angqvist2019icet,
  author  = {\\AA ngqvist, Magnus and Mu\\~{n}oz, William A. and Rahm, J. Magnus and Fransson, Erik and Durniak, C\\'eline and Rozyczko, Piotr and Rod, Thomas H. and Erhart, Paul},
  title   = {{ICET} -- A {P}ython Library for Constructing and Sampling Alloy Cluster Expansions},
  journal = {Advanced Theory and Simulations},
  year    = {2019},
  volume  = {2},
  number  = {7},
  pages   = {1900015},
  doi     = {10.1002/adts.201900015}
}""")

# The Tersoff-Hamann approximation itself (tunneling current proportional
# to the local density of states at the tip position) -- the physical
# model stb-stm's simulated STM images are built on.
TERSOFF_HAMANN = ("Tersoff1985", """@article{Tersoff1985,
  author  = {Tersoff, J. and Hamann, D. R.},
  title   = {Theory of the scanning tunneling microscope},
  journal = {Physical Review B},
  year    = {1985},
  volume  = {31},
  number  = {2},
  pages   = {805--813},
  doi     = {10.1103/PhysRevB.31.805}
}""")

# pyxtal, the random symmetry-constrained structure generation/search
# library stb-crystalcast uses for every one of its modes (generate,
# --substitute, --subgroup, --supergroup).
PYXTAL = ("Fredericks2021pyxtal", """@article{Fredericks2021pyxtal,
  author  = {Fredericks, Scott and Parrish, Kevin and Sayre, Dean and Zhu, Qiang},
  title   = {{PyXtal}: A {P}ython library for crystal structure generation and symmetry analysis},
  journal = {Computer Physics Communications},
  year    = {2021},
  volume  = {261},
  pages   = {107810},
  doi     = {10.1016/j.cpc.2020.107810}
}""")


def _read_existing_entries(bib_path):
    """Returns {key: full_entry_text} already in `bib_path`, or {} if the
    file doesn't exist -- parsed back from the same blank-line-separated
    `@type{KEY, ...}` blocks write_bib_file itself writes."""
    if not os.path.exists(bib_path):
        return {}
    with open(bib_path) as f_bib:
        content = f_bib.read()
    entries = {}
    for block in content.split("\n\n"):
        block = block.strip()
        match = re.match(r"@\w+\{([^,]+),", block)
        if match:
            entries[match.group(1)] = block
    return entries


def write_bib_file(bib_path, entries):
    """Writes a plain .bib file with one blank-line-separated entry per
    reference -- a standard BibTeX file, importable as-is into any reference
    manager (Zotero, JabRef, LaTeX's own bibtex/biber).

    Merges with whatever is already at `bib_path` (by citation key) instead
    of overwriting it -- multiple stb tools are often run in the same
    directory against complementary outputs (e.g. stb-inputfile's calc.fdf
    and stb-kpath's kpath_bs.fdf, for a 'bands' calculation), and each
    should be able to add its own citations without erasing the other's.
    Found live: stb-kpath's own write was silently dropping stb-inputfile's
    PBE/Monkhorst-Pack entries when both ran in the same output/bands/
    folder in examples/1.1-stb-inputfile/example_1.1.sh.
    """
    combined = _read_existing_entries(bib_path)
    for key, entry in entries:
        combined[key] = entry
    with open(bib_path, "w") as f_bib:
        f_bib.write("% References for the SIESTA-related calculation(s) run in this directory.\n")
        f_bib.write("% Cite these alongside the SIESTA code itself.\n\n")
        for entry in combined.values():
            f_bib.write(entry + "\n\n")
