# PseudoDojo: norm-conserving (ONCVPSP), PSML format (for SIESTA)

Downloaded from www.pseudo-dojo.org and organized on 17/11/2026.
Folder names follow `<relativity>_<xc>_<version>_<accuracy>`. The suite names each table
`dojo-<relativity>-<xc>-<accuracy>` (see `core/pseudopotentials.py`).

| Folder | `-p` name | Source table | Relativity | XC | Elements |
|---|---|---|---|---|---|
| `sr_pbe_0.5_standard` | `dojo` (= `dojo-sr-pbe-standard`) | nc-sr-05 pbe standard        | scalar | PBE        | 72 |
| `sr_pbe_0.5_stringent` | `dojo-stringent` | nc-sr-05 pbe stringent       | scalar | PBE        | 72 |
| `sr_pbesol_0.4_standard` | `dojo-pbesol` | nc-sr-04 pbesol standard     | scalar | PBEsol     | 72 |
| `sr_pbesol_0.4_stringent` | `dojo-pbesol-stringent` | nc-sr-04 pbesol stringent    | scalar | PBEsol     | 72 |
| `sr_pw_0.4_standard` | `dojo-pw` | nc-sr-04 pw standard         | scalar | LDA (PW92) | 70 |
| `sr_pw_0.4_stringent` | `dojo-pw-stringent` | nc-sr-04 pw stringent        | scalar | LDA (PW92) | 70 |
| `fr_pbe_0.4_standard` | `dojo-fr` | nc-fr-04 pbe standard        | full (SOC) | PBE    | 70 (no La or Lu) |
| `fr_pbe_0.4_stringent` | `dojo-fr-stringent` | nc-fr-04 pbe stringent       | full (SOC) | PBE    | 72 |
| `fr_pbesol_0.4_standard` | `dojo-fr-pbesol` | nc-fr-04 pbesol standard     | full (SOC) | PBEsol | 71 (no Lu) |
| `fr_pbesol_0.4_stringent` | `dojo-fr-pbesol-stringent` | nc-fr-04 pbesol stringent    | full (SOC) | PBEsol | 71 (no Lu) |
| `sr3plus_pbe_0.4_standard` | `dojo-sr3plus` | nc-sr-04-3plus pbe standard  | scalar | PBE        | 14 lanthanides |

- The folder's XC must match the calculation's `XC.Functional`/`XC.Authors`.
- **Standard or stringent:** standard is the default choice.
  - **Groups 13–16 (Ga–Se, In–Te, Tl–Po):** standard **already has the d shell in the valence**; Ge and Sn have z = 14. Stringent adds the semicore s and p of the shell below; Ge and Sn have z = 22.
  - **Hf and Ta:** stringent adds the 4f shell (z goes from 12 to 26 and from 13 to 27).
  - More electrons and more orbitals make the calculation more expensive.
  - To check a file: `grep -o 'z-pseudo="[0-9.]*"' <El>.psml`.
  - *Correction note, 18/11/2026: an earlier version of this item said stringent adds d semicore states. That was wrong.*
- **fr_*:**
  - Required for `Spin spin-orbit`. The `sr_*` files abort with `The PSML file does not contain lj projs`.
  - They also run without SOC, but do not mix absolute energies from `fr` and `sr` files.
  - There is no FR table for LDA.
- **sr3plus_*:** lanthanides Ce–Lu (no La) with the 4f electrons frozen in the core. Use them only for Ln3+ compounds. They do not describe f physics (magnetism, 4f states).
- Usage in stb: `-p dojo[-sr|-fr|-sr3plus][-pbe|-pbesol|-pw][-standard|-stringent]` (words in any order,
  defaults sr, pbe, standard); `stb-pseudo --list-banks` lists them. Every tool that copies these files into
  a run folder warns when the folder's .fdf asks for another functional or for spin-orbit with `sr` files.

Reference: M. J. van Setten et al., Comput. Phys. Commun. 226, 39 (2018), doi:10.1016/j.cpc.2018.01.012.
