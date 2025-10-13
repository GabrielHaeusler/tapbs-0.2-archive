# Third-Party Notices

This repository vendors two external open-source components used historically by TAPBS 0.2.
All third-party code is preserved with original license notices intact.

---

## APBS 1.3.0

**Name:** Adaptive Poisson–Boltzmann Solver (APBS)  
**Version:** 1.3.0  
**Upstream:** E. Baker et al., University of Washington  
**Source URL (archived):** [https://launchpad.net/ubuntu/+source/apbs/1.3.0-2](https://launchpad.net/ubuntu/+source/apbs/1.3.0-2)  
**License:** BSD-like (APBS-specific license in `LICENSE` within the source tree)  
**Retreived:** 2025-10-13
**Purpose:** Solves the linearized Poisson–Boltzmann equation (LPBE) to compute electrostatic potentials and energies.


### Patch files included in this repository

The following patches modify APBS 1.3 internals to match the workflow expected by TAPBS:

Patched files (added by TAPBS):
- `src/patch_vpbec`
- `src/patch_vpbeh`
- `src/patch_vpmgc`
- `src/patch_vpmgh`

- **`src/patch_vpbec` / `src/patch_vpbeh`**  
  Adds `Vpbe_run_ctor`/`Vpbe_run_ctor2` and extra membrane parameters; computes ionic strength, Debye length, and accessibility setup in a constructor path callable directly by TAPBS.

- **`src/patch_vpmgc` / `src/patch_vpmgh`**  
  Introduces `Vpmg_map_*` and `Vpmg_run_*` constructor/fill/destructor variants that:
  - reduce memory consumption by skipping arrays not needed for accessibility maps,  
  - allow using **precomputed dielectric/kappa/charge maps**,  
  - provide focused-boundary handling and external-energy accounting,  
  - remove verbose debugging and one redundant memory free.

These patches are applied during the APBS 1.3 build step in the README instructions.

---

## MALOC 1.5.0

**Name:** Minimal Abstraction Layer for Object-oriented C (MALOC)  
**Version:** 1.5.0  
**Upstream:** M. Holst, The FETK Project, UC San Diego  
**Source URL:** [http://www.fetk.org/codes/maloc.html](http://www.fetk.org/codes/maloc.html)  
**License:** FETK/MALOC License (BSD-like)  
**Retreived:** 2025-10-13
**Purpose:** Foundation library used by APBS for memory management, timing, and portable abstraction across architectures.

---

## Notes

Both APBS and MALOC were obtained from publicly available sources that permit redistribution.  
They are included here **solely for archival reproducibility**, under their respective open-source licenses.  
All original copyright headers and license files remain unmodified within their directories.

If you are a rightsholder and wish to clarify licensing or request changes, please open an issue in this repository.