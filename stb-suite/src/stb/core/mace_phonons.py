"""Shared MACE-MP-0 finite-difference phonon force-constant logic, used by
stb-phononsML and stb-phononsQHA. Callers must call core.deps.require_mace()
themselves first, same convention as core/mace_relax.py -- this module only
imports ase/phonopy/mace lazily inside the functions below, never at module
level.
"""

import numpy as np


def generate_ml_displacements(atoms, dim, distance, calc, relax=True, fmax=0.05,
                               max_steps=200):
    """Builds a Phonopy object for `atoms` (an ase.Atoms, positions/cell
    already in Angstrom) and generates its finite-difference displacement
    dataset -- the MACE-only equivalent of phonons_create.py's structure
    read + generate_displacements step, minus any SIESTA/bohr involvement.

    Builds the Phonopy object with calculator=None (the default Angstrom/
    eV/eV-per-Angstrom convention, the same one ASE/MACE already use) --
    unlike this workflow's SIESTA path, which keeps structures internally
    in bohr (see phonons_create.py), no unit conversion is needed here.

    `relax`, if True (recommended -- a phonon calculation assumes ~zero net
    force at the reference geometry), relaxes ATOMIC POSITIONS ONLY before
    generating displacements; the cell is left exactly as given, since
    callers may have deliberately fixed it (e.g. stb-phononsQHA locks the
    cell to a specific target volume per point -- relaxing it here would
    defeat that).

    Returns (phonon, atoms, relax_info): `phonon.supercells_with_displacements`
    is ready, but no forces/force constants yet -- see compute_force_constants
    below. `atoms` is the (possibly relaxed) ase.Atoms actually used to build
    the Phonopy unit cell, calculator still attached -- callers that also
    need the electronic energy at this geometry (e.g. stb-phononsQHA's
    per-volume total energy) can call atoms.get_potential_energy() on it
    without a redundant MACE evaluation. relax_info is None when relax=False,
    otherwise (f0, f1, converged, steps_used) -- residual force before/after,
    optimizer convergence flag, and step count, for the caller to report.
    """
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from stb.core import mace_relax

    atoms = atoms.copy()
    atoms.calc = calc

    relax_info = None
    if relax:
        f0 = np.abs(atoms.get_forces()).max()
        converged, steps_used = mace_relax.relax(
            atoms, calc, cell_mask=None, fmax=fmax, max_steps=max_steps)
        f1 = np.abs(atoms.get_forces()).max()
        relax_info = (f0, f1, converged, steps_used)

    unitcell = PhonopyAtoms(numbers=atoms.numbers, positions=atoms.get_positions(),
                             cell=atoms.get_cell())
    supercell_matrix = [[dim[0], 0, 0], [0, dim[1], 0], [0, 0, dim[2]]]
    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(distance=distance)

    return phonon, atoms, relax_info


def compute_force_constants(phonon, calc, progress_callback=None):
    """Computes MACE forces for every supercell in
    `phonon.supercells_with_displacements` and builds force constants from
    them, in place (phonon.forces set, phonon.produce_force_constants()
    called). `phonon` must already have a displacement dataset (see
    generate_ml_displacements above).

    `progress_callback(i, n)`, if given, is called after each of the `n`
    displaced supercells' forces are computed (i is 1-based) -- callers
    print their own progress line, this module stays print-free (same
    convention as core/mace_relax.py).
    """
    from ase import Atoms

    supercells = phonon.supercells_with_displacements
    force_sets = []
    for i, scell in enumerate(supercells):
        disp_atoms = Atoms(numbers=scell.numbers, positions=scell.positions,
                            cell=scell.cell, pbc=True)
        disp_atoms.calc = calc
        force_sets.append(disp_atoms.get_forces())
        if progress_callback is not None:
            progress_callback(i + 1, len(supercells))

    phonon.forces = force_sets
    phonon.produce_force_constants()
