# PseudoDojo: norm-conserving (ONCVPSP), PSML format (for SIESTA)

Downloaded from www.pseudo-dojo.org (the `<table>_psml.tgz` archives) and organized on 17/11/2026;
the two PBE tables were replaced by nc-sr-04 on 20/11/2026 (see below). Every folder has a
`manifest.json`: source archive and its SHA-256, and per file the md5, valence (`z-pseudo`),
relativity, and whether it has lj and scalar-relativistic projectors. On 20/11/2026 every file was
checked byte by byte against its official archive. `stb-pseudo --list-banks` re-checks the md5s.
Folder names follow `<relativity>_<xc>_<version>_<accuracy>`. The suite names each table
`dojo-<relativity>-<xc>-<accuracy>` (see `core/pseudopotentials.py`).

| Folder | `-p` name | Source table | Relativity | XC | Elements |
|---|---|---|---|---|---|
| `sr_pbe_0.4_standard` | `dojo` (= `dojo-sr-pbe-standard`) | nc-sr-04 pbe standard        | scalar | PBE        | 72 |
| `sr_pbe_0.4_stringent` | `dojo-stringent` | nc-sr-04 pbe stringent       | scalar | PBE        | 72 |
| `sr_pbesol_0.4_standard` | `dojo-pbesol` | nc-sr-04 pbesol standard     | scalar | PBEsol     | 72 |
| `sr_pbesol_0.4_stringent` | `dojo-pbesol-stringent` | nc-sr-04 pbesol stringent    | scalar | PBEsol     | 72 |
| `sr_pw_0.4_standard` | `dojo-pw` | nc-sr-04 pw standard         | scalar | LDA (PW92) | 70 |
| `sr_pw_0.4_stringent` | `dojo-pw-stringent` | nc-sr-04 pw stringent        | scalar | LDA (PW92) | 70 |
| `fr_pbe_0.4_standard` | `dojo-fr` | nc-fr-04 pbe standard        | full (SOC) | PBE    | 70 (no La or Lu) |
| `fr_pbe_0.4_stringent` | `dojo-fr-stringent` | nc-fr-04 pbe stringent       | full (SOC) | PBE    | 72 |
| `fr_pbesol_0.4_standard` | `dojo-fr-pbesol` | nc-fr-04 pbesol standard     | full (SOC) | PBEsol | 71 (no Lu) |
| `fr_pbesol_0.4_stringent` | `dojo-fr-pbesol-stringent` | nc-fr-04 pbesol stringent    | full (SOC) | PBEsol | 71 (no Lu) |
| `sr3plus_pbe_0.4_standard` | `dojo-sr3plus` | nc-sr-04-3plus pbe standard  | scalar | PBE        | 14 lanthanides |

- **Why nc-sr-04 and not nc-sr-05 for PBE.** Up to 2.0.0, `dojo` and `dojo-stringent` were
  nc-sr-05, which is nc-sr-04 with 11 files regenerated: Ba, Bi, I, Pb, Po, Rb, Rn, S, Te, Tl, Xe.
  Those 11 have KB projectors with `ekb="0"`. With SIESTA 5.4.2 they
  - give wrong equilibrium structures with any number of MPI ranks, with no warning (Pb fcc
    standard 14.6 against 30.8 Ang^3/atom, S2 1.65 against 1.92 Ang);
  - abort in parallel with `propor: ERROR: IMAX = 0` (from 4 ranks for some elements).

  The other 61 files of each table are byte-identical in both versions, so only these 11 changed.
  The change moves the valence of three elements:

  | Element | standard nc-sr-05 -> nc-sr-04 | stringent nc-sr-05 -> nc-sr-04 |
  |---|---|---|
  | I  | 17 -> 7 | 25 -> 17 |
  | Rn | 18 -> 18 | 26 -> 18 |
  | Xe | 8 -> 8 | 26 -> 8 |

  Results for I (both tables) and stringent Rn and Xe are therefore not comparable with runs made
  with 2.0.0. Every tool refuses a `.psml` with `ekb="0"` projectors, also one from a user folder
  (`-p <folder>`): it stops with an `[ERROR]` before the file is copied.
- The folder's XC must match the calculation's `XC.Functional`/`XC.Authors`.
- **Standard or stringent:** standard is the default choice.
  - **Groups 13–16 (Ga–Se, In–Te, Tl–Po):** standard **already has the d shell in the valence**; Ge and Sn have z = 14. Stringent adds the semicore s and p of the shell below; Ge and Sn have z = 22.
  - **Hf and Ta:** stringent adds the 4f shell (z goes from 12 to 26 and from 13 to 27).
  - More electrons and more orbitals make the calculation more expensive.
  - To check a file: `grep -o 'z-pseudo="[0-9.]*"' <El>.psml`.
  - *Correction note, 18/11/2026: an earlier version of this item said stringent adds d semicore states. That was wrong.*
- **fr_*:**
  - Required for `Spin spin-orbit`. The `sr_*` files abort with `The PSML file does not contain lj projs`.
  - Without SOC: the `fr_pbe` files run (do not mix absolute energies from `fr` and `sr` files).
    Of the `fr_pbesol` files only Ge, Sb and Te do; the other 68 have no scalar-relativistic
    projector set, and SIESTA 5.4.2 stops with 'generation of scalar-relativistic part is not
    available yet'. The tools refuse that combination with an `[ERROR]`; use `-p dojo-pbesol`.
  - There is no FR table for LDA.
- **sr3plus_*:** lanthanides Ce–Lu (no La) with the 4f electrons frozen in the core. Use them only for Ln3+ compounds. They do not describe f physics (magnetism, 4f states).
- Usage in stb: `-p dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent]` (words in any order,
  defaults sr, pbe, standard); `stb-pseudo --list-banks` lists them. Every tool that copies these files into
  a run folder warns when the folder's .fdf asks for another functional or for spin-orbit with `sr` files.

Reference: M. J. van Setten et al., Comput. Phys. Commun. 226, 39 (2018), doi:10.1016/j.cpc.2018.01.012.
