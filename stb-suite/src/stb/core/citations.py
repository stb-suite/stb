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

# The generalized Inverse Participation Ratio (IPR) paper -- the same
# citation sisl's own EigenstateElectron.ipr() docstring points to
# (":cite:`Murphy2011`") for the Anderson-localization criterion stb-ipr
# reports.
IPR_MURPHY = ("Murphy2011", """@article{Murphy2011,
  author  = {Murphy, N. C. and Wortis, R. and Atkinson, W. A.},
  title   = {Generalized inverse participation ratio as a possible measure of localization for interacting systems},
  journal = {Physical Review B},
  year    = {2011},
  volume  = {83},
  number  = {18},
  pages   = {184206},
  doi     = {10.1103/PhysRevB.83.184206}
}""")

# The original Crystal Orbital Overlap Population (COOP) paper -- the
# bonding/antibonding curve stb-coop's --quantity coop reproduces.
COOP_HUGHBANKS = ("Hughbanks1983coop", """@article{Hughbanks1983coop,
  author  = {Hughbanks, Timothy and Hoffmann, Roald},
  title   = {Chains of trans-edge-sharing molybdenum octahedra: metal-metal bonding in extended systems},
  journal = {Journal of the American Chemical Society},
  year    = {1983},
  volume  = {105},
  number  = {11},
  pages   = {3528--3537},
  doi     = {10.1021/ja00349a027}
}""")

# The Crystal Orbital Hamilton Population (COHP) paper -- stb-coop's
# --quantity cohp.
COHP_DRONSKOWSKI = ("Dronskowski1993cohp", """@article{Dronskowski1993cohp,
  author  = {Dronskowski, Richard and Bl\\"ochl, Peter E.},
  title   = {Crystal Orbital Hamilton Populations ({COHP}): Energy-Resolved Visualization of Chemical Bonding in Solids Based on Density-Functional Calculations},
  journal = {The Journal of Physical Chemistry},
  year    = {1993},
  volume  = {97},
  number  = {33},
  pages   = {8617--8624},
  doi     = {10.1021/j100135a014}
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


# Pseudopotential-bank citations, one per stb.core.pseudopotentials.BANKS
# key -- kept here rather than read from BANKS[...]['citation'] directly
# since that field is a plain human-readable string (printed by
# resolve_pseudo_source()), not already-formatted BibTeX. Moved here from
# inputfile.py once stb-pseudo became a second consumer.
PSEUDO_DOJO = ("vanSetten2018", """@article{vanSetten2018,
  author  = {van Setten, M. J. and Giantomassi, M. and Bousquet, E. and Verstraete, M. J. and Hamann, D. R. and Gonze, X. and Rignanese, G.-M.},
  title   = {The {PseudoDojo}: Training and grading a 85 element optimized norm-conserving pseudopotential table},
  journal = {Computer Physics Communications},
  year    = {2018},
  volume  = {226},
  pages   = {39--54},
  doi     = {10.1016/j.cpc.2018.01.012}
}""")

PSEUDO_VIRTUAL_VAULT = ("VirtualVault", """@misc{VirtualVault,
  author       = {{NNIN/C}},
  title        = {{SIESTA} Pseudopotentials Virtual Vault},
  institution  = {Cornell NanoScale Science and Technology Facility},
  url          = {https://nninc.cnf.cornell.edu/}
}""")

# Keyed by stb.core.pseudopotentials.BANKS name, for callers that already
# have the resolved bank name (e.g. via detect_pseudo_bank()) and just want
# the matching bib entry.
PSEUDO_BANK_CITATIONS = {
    "dojo": PSEUDO_DOJO,
    "virtual_vault": PSEUDO_VIRTUAL_VAULT,
}



# ---- Entries copied verbatim from the verified bibliography of the STB-SUITE paper. ----
PYMATGEN = ("Ong2013pymatgen", """@article{Ong2013pymatgen,
  author  = {Ong, Shyue Ping and Richards, William Davidson and Jain, Anubhav and others},
  title   = {{Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis}},
  journal = {Computational Materials Science},
  year    = {2013},
  volume  = {68},
  pages   = {314--319},
  doi     = {10.1016/j.commatsci.2012.10.028}
}""")
ASE = ("Larsen2017ase", """@article{Larsen2017ase,
  author  = {Larsen, Ask Hjorth and Mortensen, Jens J{\\o}rgen and Blomqvist, Jakob and others},
  title   = {{The atomic simulation environment---a Python library for working with atoms}},
  journal = {Journal of Physics: Condensed Matter},
  year    = {2017},
  volume  = {29},
  number  = {27},
  pages   = {273002},
  doi     = {10.1088/1361-648X/aa680e}
}""")
PHONOPY = ("Togo2015phonopy", """@article{Togo2015phonopy,
  author  = {Togo, Atsushi and Tanaka, Isao},
  title   = {{First principles phonon calculations in materials science}},
  journal = {Scripta Materialia},
  year    = {2015},
  volume  = {108},
  pages   = {1--5},
  doi     = {10.1016/j.scriptamat.2015.07.021}
}""")
PHONOPY_IMPL = ("Togo2023phonopy", """@article{Togo2023phonopy,
  author  = {Togo, Atsushi and Chaput, Laurent and Tadano, Terumasa and Tanaka, Isao},
  title   = {{Implementation strategies in phonopy and phono3py}},
  journal = {Journal of Physics: Condensed Matter},
  year    = {2023},
  volume  = {35},
  number  = {35},
  pages   = {353001},
  doi     = {10.1088/1361-648X/acd831}
}""")
SPGLIB = ("Togo2024spglib", """@article{Togo2024spglib,
  author  = {Togo, Atsushi and Shinohara, Kohei and Tanaka, Isao},
  title   = {{Spglib: a software library for crystal symmetry search}},
  journal = {Science and Technology of Advanced Materials: Methods},
  year    = {2024},
  volume  = {4},
  number  = {1},
  pages   = {2384822},
  doi     = {10.1080/27660400.2024.2384822}
}""")
SEEKPATH = ("Hinuma2017seekpath", """@article{Hinuma2017seekpath,
  author  = {Hinuma, Yoyo and Pizzi, Giovanni and Kumagai, Yu and Oba, Fumiyasu and Tanaka, Isao},
  title   = {{Band structure diagram paths based on crystallography}},
  journal = {Computational Materials Science},
  year    = {2017},
  volume  = {128},
  pages   = {140--184},
  doi     = {10.1016/j.commatsci.2016.10.015}
}""")
SISL = ("Papior2017sisl", """@misc{Papior2017sisl,
  author    = {Papior, Nick and Febrer Calabozo, Pol},
  title     = {sisl},
  year      = {2026},
  publisher = {Zenodo},
  note      = {Software, concept DOI (year of the latest release); cite the release used},
  doi       = {10.5281/zenodo.597181}
}""")
BADER_HENKELMAN = ("Henkelman2006bader", """@article{Henkelman2006bader,
  author  = {Henkelman, Graeme and Arnaldsson, Andri and J{\\'o}nsson, Hannes},
  title   = {{A fast and robust algorithm for Bader decomposition of charge density}},
  journal = {Computational Materials Science},
  year    = {2006},
  volume  = {36},
  number  = {3},
  pages   = {354--360},
  doi     = {10.1016/j.commatsci.2005.04.010}
}""")
BADER_TANG = ("Tang2009bader", """@article{Tang2009bader,
  author  = {Tang, W and Sanville, E and Henkelman, G},
  title   = {{A grid-based Bader analysis algorithm without lattice bias}},
  journal = {Journal of Physics: Condensed Matter},
  year    = {2009},
  volume  = {21},
  number  = {8},
  pages   = {084204},
  doi     = {10.1088/0953-8984/21/8/084204}
}""")
CINEB = ("Henkelman2000cineb", """@article{Henkelman2000cineb,
  author  = {Henkelman, Graeme and Uberuaga, Blas P. and J{\\'o}nsson, Hannes},
  title   = {{A climbing image nudged elastic band method for finding saddle points and minimum energy paths}},
  journal = {The Journal of Chemical Physics},
  year    = {2000},
  volume  = {113},
  number  = {22},
  pages   = {9901--9904},
  doi     = {10.1063/1.1329672}
}""")
IDPP = ("Smidstrup2014idpp", """@article{Smidstrup2014idpp,
  author  = {Smidstrup, S{\\o}ren and Pedersen, Andreas and Stokbro, Kurt and J{\\'o}nsson, Hannes},
  title   = {{Improved initial guess for minimum energy path calculations}},
  journal = {The Journal of Chemical Physics},
  year    = {2014},
  volume  = {140},
  number  = {21},
  pages   = {214106},
  doi     = {10.1063/1.4878664}
}""")
DFT_D3 = ("Grimme2010d3", """@article{Grimme2010d3,
  author  = {Grimme, Stefan and Antony, Jens and Ehrlich, Stephan and Krieg, Helge},
  title   = {{A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu}},
  journal = {The Journal of Chemical Physics},
  year    = {2010},
  volume  = {132},
  number  = {15},
  pages   = {154104},
  doi     = {10.1063/1.3382344}
}""")
BSSE_BOYS = ("Boys1970bsse", """@article{Boys1970bsse,
  author  = {Boys, S.F. and Bernardi, F.},
  title   = {{The calculation of small molecular interactions by the differences of separate total energies. Some procedures with reduced errors}},
  journal = {Molecular Physics},
  year    = {1970},
  volume  = {19},
  number  = {4},
  pages   = {553--566},
  doi     = {10.1080/00268977000101561}
}""")
BIRCH = ("Birch1947eos", """@article{Birch1947eos,
  author  = {Birch, Francis},
  title   = {{Finite Elastic Strain of Cubic Crystals}},
  journal = {Physical Review},
  year    = {1947},
  volume  = {71},
  number  = {11},
  pages   = {809--824},
  doi     = {10.1103/physrev.71.809}
}""")
HIRSHFELD = ("Hirshfeld1977", """@article{Hirshfeld1977,
  author  = {Hirshfeld, F. L.},
  title   = {{Bonded-atom fragments for describing molecular charge densities}},
  journal = {Theoretica Chimica Acta},
  year    = {1977},
  volume  = {44},
  number  = {2},
  pages   = {129--138},
  doi     = {10.1007/bf00549096}
}""")
HIRSHFELD_I = ("Bultinck2007hirshfeldI", """@article{Bultinck2007hirshfeldI,
  author  = {Bultinck, Patrick and Van Alsenoy, Christian and Ayers, Paul W. and Carb{\\'o}-Dorca, Ramon},
  title   = {{Critical analysis and extension of the Hirshfeld atoms in molecules}},
  journal = {The Journal of Chemical Physics},
  year    = {2007},
  volume  = {126},
  number  = {14},
  pages   = {144111},
  doi     = {10.1063/1.2715563}
}""")
HUBBARD_LR = ("Cococcioni2005hubbard", """@article{Cococcioni2005hubbard,
  author  = {Cococcioni, Matteo and de Gironcoli, Stefano},
  title   = {{Linear response approach to the calculation of the effective interaction parameters in the LDA + U method}},
  journal = {Physical Review B},
  year    = {2005},
  volume  = {71},
  number  = {3},
  pages   = {035105},
  doi     = {10.1103/PhysRevB.71.035105}
}""")
SQS = ("Zunger1990sqs", """@article{Zunger1990sqs,
  author  = {Zunger, Alex and Wei, S.-H. and Ferreira, L. G. and Bernard, James E.},
  title   = {{Special quasirandom structures}},
  journal = {Physical Review Letters},
  year    = {1990},
  volume  = {65},
  number  = {3},
  pages   = {353--356},
  doi     = {10.1103/PhysRevLett.65.353}
}""")
ZSL = ("Zur1984zsl", """@article{Zur1984zsl,
  author  = {Zur, A. and McGill, T. C.},
  title   = {{Lattice match: An application to heteroepitaxy}},
  journal = {Journal of Applied Physics},
  year    = {1984},
  volume  = {55},
  number  = {2},
  pages   = {378--386},
  doi     = {10.1063/1.333084}
}""")
SETYAWAN = ("Setyawan2010kpath", """@article{Setyawan2010kpath,
  author  = {Setyawan, Wahyu and Curtarolo, Stefano},
  title   = {{High-throughput electronic band structure calculations: Challenges and tools}},
  journal = {Computational Materials Science},
  year    = {2010},
  volume  = {49},
  number  = {2},
  pages   = {299--312},
  doi     = {10.1016/j.commatsci.2010.05.010}
}""")
MATPROJ = ("Jain2013matproj", """@article{Jain2013matproj,
  author  = {Jain, Anubhav and Ong, Shyue Ping and Hautier, Geoffroy and others},
  title   = {{Commentary: The Materials Project: A materials genome approach to accelerating materials innovation}},
  journal = {APL Materials},
  year    = {2013},
  volume  = {1},
  number  = {1},
  pages   = {011002},
  doi     = {10.1063/1.4812323}
}""")
OPTIMADE = ("Andersen2021optimade", """@article{Andersen2021optimade,
  author  = {Andersen, Casper W. and Armiento, Rickard and Blokhin, Evgeny and others},
  title   = {{OPTIMADE, an API for exchanging materials data}},
  journal = {Scientific Data},
  year    = {2021},
  volume  = {8},
  number  = {1},
  pages   = {217},
  doi     = {10.1038/s41597-021-00974-z}
}""")
COD = ("Grazulis2012cod", """@article{Grazulis2012cod,
  author  = {Gra{\\v z}ulis, Saulius and Da{\\v s}kevi{\\v c}, Adriana and Merkys, Andrius and others},
  title   = {{Crystallography Open Database (COD): an open-access collection of crystal structures and platform for world-wide collaboration}},
  journal = {Nucleic Acids Research},
  year    = {2011},
  volume  = {40},
  number  = {D1},
  pages   = {D420--D427},
  doi     = {10.1093/nar/gkr900}
}""")
PSML = ("Garcia2018psml", """@article{Garcia2018psml,
  author  = {Garc{\\'\\i}a, Alberto and Verstraete, Matthieu J. and Pouillon, Yann and Junquera, Javier},
  title   = {{The psml format and library for norm-conserving pseudopotential data curation and interoperability}},
  journal = {Computer Physics Communications},
  year    = {2018},
  volume  = {227},
  pages   = {51--71},
  doi     = {10.1016/j.cpc.2018.02.011}
}""")
WANTIBEXOS = ("Dias2023wantibexos", """@article{Dias2023wantibexos,
  author  = {Dias, Alexandre C. and Silveira, Julian F.R.V. and Qu, Fanyao},
  title   = {{WanTiBEXOS: A Wannier based Tight Binding code for electronic band structure, excitonic and optoelectronic properties of solids}},
  journal = {Computer Physics Communications},
  year    = {2023},
  volume  = {285},
  pages   = {108636},
  doi     = {10.1016/j.cpc.2022.108636}
}""")

# ---- Further standard references (library and method papers). ----
NUMPY = ("Harris2020numpy", """@article{Harris2020numpy,
  author  = {Harris, Charles R. and Millman, K. Jarrod and van der Walt, St{\\'e}fan J. and others},
  title   = {{Array programming with NumPy}},
  journal = {Nature},
  year    = {2020},
  volume  = {585},
  pages   = {357--362},
  doi     = {10.1038/s41586-020-2649-2}
}""")
SCIPY = ("Virtanen2020scipy", """@article{Virtanen2020scipy,
  author  = {Virtanen, Pauli and Gommers, Ralf and Oliphant, Travis E. and others},
  title   = {{SciPy 1.0: fundamental algorithms for scientific computing in Python}},
  journal = {Nature Methods},
  year    = {2020},
  volume  = {17},
  pages   = {261--272},
  doi     = {10.1038/s41592-019-0686-2}
}""")
MATPLOTLIB = ("Hunter2007matplotlib", """@article{Hunter2007matplotlib,
  author  = {Hunter, John D.},
  title   = {{Matplotlib: A 2D graphics environment}},
  journal = {Computing in Science \\& Engineering},
  year    = {2007},
  volume  = {9},
  number  = {3},
  pages   = {90--95},
  doi     = {10.1109/MCSE.2007.55}
}""")
PYTORCH = ("Paszke2019pytorch", """@inproceedings{Paszke2019pytorch,
  author    = {Paszke, Adam and Gross, Sam and Massa, Francisco and others},
  title     = {{PyTorch: An Imperative Style, High-Performance Deep Learning Library}},
  booktitle = {Advances in Neural Information Processing Systems 32 (NeurIPS 2019)},
  year      = {2019},
  pages     = {8024--8035},
  eprint    = {1912.01703},
  archivePrefix = {arXiv}
}""")
GQCA_SHER = ("Sher1987gqca", """@article{Sher1987gqca,
  author  = {Sher, A. and van Schilfgaarde, M. and Chen, A.-B. and Chen, W.},
  title   = {{Quasichemical approximation in binary alloys}},
  journal = {Physical Review B},
  year    = {1987},
  volume  = {36},
  number  = {8},
  pages   = {4279--4295},
  doi     = {10.1103/PhysRevB.36.4279}
}""")
MURNAGHAN = ("Murnaghan1944", """@article{Murnaghan1944,
  author  = {Murnaghan, F. D.},
  title   = {{The Compressibility of Media under Extreme Pressures}},
  journal = {Proceedings of the National Academy of Sciences},
  year    = {1944},
  volume  = {30},
  number  = {9},
  pages   = {244--247},
  doi     = {10.1073/pnas.30.9.244}
}""")
POREZAG = ("Porezag1996", """@article{Porezag1996,
  author  = {Porezag, Dirk and Pederson, Mark R.},
  title   = {{Infrared intensities and Raman-scattering activities within density-functional theory}},
  journal = {Physical Review B},
  year    = {1996},
  volume  = {54},
  number  = {11},
  pages   = {7830--7836},
  doi     = {10.1103/PhysRevB.54.7830}
}""")
LEPAGE_SAXE = ("LePage2002", """@article{LePage2002,
  author  = {Le Page, Yvon and Saxe, Paul},
  title   = {{Symmetry-general least-squares extraction of elastic data for strained materials from ab initio calculations of stress}},
  journal = {Physical Review B},
  year    = {2002},
  volume  = {65},
  number  = {10},
  pages   = {104104},
  doi     = {10.1103/PhysRevB.65.104104}
}""")
HILL = ("Hill1952", """@article{Hill1952,
  author  = {Hill, R.},
  title   = {{The Elastic Behaviour of a Crystalline Aggregate}},
  journal = {Proceedings of the Physical Society. Section A},
  year    = {1952},
  volume  = {65},
  number  = {5},
  pages   = {349--354},
  doi     = {10.1088/0370-1298/65/5/307}
}""")
MOUHAT_COUDERT = ("Mouhat2014", """@article{Mouhat2014,
  author  = {Mouhat, F{\\'e}lix and Coudert, Fran{\\c{c}}ois-Xavier},
  title   = {{Necessary and sufficient elastic stability conditions in various crystal systems}},
  journal = {Physical Review B},
  year    = {2014},
  volume  = {90},
  number  = {22},
  pages   = {224104},
  doi     = {10.1103/PhysRevB.90.224104}
}""")
DIPOLE_BENGTSSON = ("Bengtsson1999", """@article{Bengtsson1999,
  author  = {Bengtsson, Lennart},
  title   = {{Dipole correction for surface supercell calculations}},
  journal = {Physical Review B},
  year    = {1999},
  volume  = {59},
  number  = {19},
  pages   = {12301--12304},
  doi     = {10.1103/PhysRevB.59.12301}
}""")

def pseudo_bank_citation(bank):
    """Bib entry of a bundled bank name -- every PseudoDojo table (dojo-sr-pbe-standard,
    dojo-fr-pbesol-stringent, ...) cites the same paper. None for an unknown name."""
    if not bank:
        return None
    if bank in PSEUDO_BANK_CITATIONS:
        return PSEUDO_BANK_CITATIONS[bank]
    if str(bank).lower().startswith("dojo"):
        return PSEUDO_DOJO
    return None


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


# ---- One references file per run, always written (suite rule) -----------------------------
#
# Every command writes, next to its report and with the same number, the BibTeX file of the
# run: stb_<tool>_references.bib in the folder it runs from (stb_<tool>_references_1.bib, ...
# when the report is stb_<tool>_report_1.txt), never overwritten. It holds
#   - what the tool cites for its method (TOOL_CITATIONS below),
#   - what the run itself cited (cite(): an XC functional, a pseudopotential bank, D3, ...),
#   - every library the run actually loaded (LIBRARY_CITATIONS, read from sys.modules at exit).
# core.cli.new_report_path registers the run; the file is written when the process exits.

# Python module (top-level name) -> entries to cite when the run loaded it.
LIBRARY_CITATIONS = {
    "numpy": [NUMPY], "scipy": [SCIPY], "matplotlib": [MATPLOTLIB],
    "pymatgen": [PYMATGEN], "ase": [ASE], "spglib": [SPGLIB], "seekpath": [SEEKPATH],
    "sisl": [SISL], "phonopy": [PHONOPY, PHONOPY_IMPL], "pybader": [BADER_HENKELMAN, BADER_TANG],
    "icet": [ICET], "pyxtal": [PYXTAL], "mace": [MACE], "torch": [PYTORCH],
}

# Methods each command implements or relies on, cited on every run of it (its libraries come
# from LIBRARY_CITATIONS; run-dependent ones -- functional, bank, D3 when forced -- via cite()).
_S = [SIESTA, SIESTA_RECENT]
TOOL_CITATIONS = {
    # 1 Inputs
    "inputfile": _S, "pseudo": _S, "kgrid": _S + [MONKHORST_PACK], "kpath": _S + [SETYAWAN],
    "dftu": _S, "fetch": [], "mlrelax": [],
    # 2 Analysis
    "bands": _S, "dos": _S, "convdos": _S, "structural": [], "symmetry": [], "bader": _S + [BADER_HENKELMAN, BADER_TANG],
    "workfunction": _S + [DIPOLE_BENGTSSON], "density": _S, "xrd": [], "fatbands": _S, "stm": _S + [TERSOFF_HAMANN],
    "wfdensity": _S, "sts": _S + [TERSOFF_HAMANN], "coop": _S + [COOP_HUGHBANKS, COHP_DRONSKOWSKI],
    "ipr": _S + [IPR_MURPHY], "effmass": _S, "spintexture": _S, "aimdAnalysis": _S, "nativecharges": _S + [HIRSHFELD],
    # 3 Utils
    "translate": [], "clean": _S, "cube": _S, "siesta2wtb": _S + [WANTIBEXOS], "ani2traj": _S, "status": _S,
    "archive": _S,
    # 4 Structures
    "2Dstacking": [ZSL], "supercell": [], "slab": [], "nanotube": [], "defect": [], "sqs": [SQS], "unitcell": [],
    "crystalbuilder": [], "passivate": [], "molecule": [], "amorphize": [], "crystalcast": [],
    # 5 Workflow
    "strain": _S, "strainAnalysis": _S,
    "elasticInputs": _S + [LEPAGE_SAXE], "elasticAnalysis": _S + [LEPAGE_SAXE, HILL, MOUHAT_COUDERT],
    "cohesive": _S + [BSSE_BOYS], "cohesiveAnalysis": _S + [BSSE_BOYS],
    "phononsCreate": _S, "phononsPos": _S,
    "convergence": _S, "convergenceAnalysis": _S,
    "xrdsearch": [], "xrdrank": [],
    "hubbardu": _S + [HUBBARD_LR], "hubbarduAlphas": _S + [HUBBARD_LR], "hubbarduAnalysis": _S + [HUBBARD_LR],
    "adsorb": _S + [DFT_D3, DIPOLE_BENGTSSON], "adsorbBsse": _S + [BSSE_BOYS], "adsorbAnalysis": _S + [BSSE_BOYS],
    "adsorbGibbs": _S,
    "neb": _S + [CINEB, IDPP], "nebCycle": _S + [CINEB], "nebAnalysis": _S + [CINEB],
    "stackingfault": _S, "stackingfaultAnalysis": _S,
    "raman": _S + [POREZAG], "ramanModes": _S + [POREZAG], "ramanAnalysis": _S + [POREZAG],
    "ir": _S + [POREZAG], "irModes": _S + [POREZAG], "irAnalysis": _S + [POREZAG],
    "gqca": _S + [GQCA_SHER], "gqcaAnalysis": _S + [GQCA_SHER],
    "optical": _S, "opticalAnalysis": _S,
    "eosInputs": _S + [BIRCH, MURNAGHAN], "eosAnalysis": _S + [BIRCH, MURNAGHAN],
    "chargediffPrep": _S + [BSSE_BOYS], "chargediffAnalysis": _S + [BSSE_BOYS],
    "hirshfeldPrep": _S + [HIRSHFELD, HIRSHFELD_I], "hirshfeldIons": _S + [HIRSHFELD, HIRSHFELD_I],
    "hirshfeldAnalysis": _S + [HIRSHFELD, HIRSHFELD_I],
    # 7 ML Simulations
    "mlff": _S + [MACE_MP], "mlffAnalysis": _S + [MACE_MP], "mlffValidate": _S, "mlffActiveLearning": _S,
    "mlmd": [], "mlphonons": [], "mlelastic": [LEPAGE_SAXE, HILL, MOUHAT_COUDERT], "mlneb": [CINEB, IDPP],
    "mladsorb": [],
}

_RUN = {"tool": None, "path": None, "entries": {}, "written": False, "announced": False}


def references_name(tool, n=0):
    """stb_<tool>_references.bib (n=0) or stb_<tool>_references_<n>.bib."""
    tool = tool[4:] if tool.startswith("stb-") else tool
    return f"stb_{tool}_references.bib" if n == 0 else f"stb_{tool}_references_{n}.bib"


def register_run(tool, path):
    """Called by core.cli.new_report_path: the references file of this run and its tool. The
    file is written when the process exits (also after an error, with what the run used)."""
    import atexit
    first = _RUN["path"] is None
    _RUN.update(tool=tool[4:] if tool.startswith("stb-") else tool, path=path)
    if first:
        atexit.register(write_run_references)


def cite(entries):
    """Adds (key, bibtex) entries this run uses to its references file. Accepts one entry or a
    list. Returns the path of the run's references file (None before the run is registered)."""
    if entries and isinstance(entries, tuple) and isinstance(entries[0], str):
        entries = [entries]
    for key, text in entries or []:
        _RUN["entries"].setdefault(key, text)
    return _RUN["path"]


def run_references_path():
    """Path of this run's references file (stb_<tool>_references[_N].bib), or None. A tool that
    shows it in its report need not be told again at exit."""
    _RUN["announced"] = True
    return _RUN["path"]


def loaded_library_entries():
    """Entries of the libraries in LIBRARY_CITATIONS that this process has imported."""
    import sys
    out = []
    for module, entries in LIBRARY_CITATIONS.items():
        if module in sys.modules:
            out.extend(entries)
    return out


def write_run_references():
    """Writes the run's references file: the tool's method entries, the run's own cite()s and
    the loaded libraries, without duplicates (by key), in that order."""
    if _RUN["path"] is None or _RUN["written"]:
        return
    tool = _RUN["tool"]
    ordered = {}
    for key, text in list(TOOL_CITATIONS.get(tool, [])) + list(_RUN["entries"].items()) + loaded_library_entries():
        ordered.setdefault(key, text)
    try:
        with open(_RUN["path"], "w") as f_bib:
            f_bib.write(f"% References of this run of stb-{tool} (STB-SUITE): the methods it uses and the\n")
            f_bib.write("% libraries it loaded. One file per run, never overwritten; cite them in your work.\n\n")
            for text in ordered.values():
                f_bib.write(text + "\n\n")
        _RUN["written"] = True
        if not _RUN["announced"]:
            print(f"References: {_RUN['path']}")
    except OSError:
        pass
