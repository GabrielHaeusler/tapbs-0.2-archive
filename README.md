# TAPBS 0.2

Front end to **APBS 1.3** that automates the many LPBE (linearized Poisson–Boltzmann) calculations needed for continuum-electrostatics **pKₐ** and **redox** workflows.
Prepares inputs and aggregates electrostatic energy terms consumed by **Karlsberg 2.0** (a.k.a. Karlsberg+).

> TAPBS supports local conformational changes, multiple charge states, regions with different dielectric constants, and implicit membrane slabs.  
> See the original HTML manual in `doc/userguide.html`.

---

## Archival scope (what this repository guarantees)

This is an **archival re-release** that reproduces the historical software stack used by TAPBS 0.2.  
It vendors the exact third-party sources required to build TAPBS **as it worked then**:

- `external/apbs-1.3-source/` — APBS 1.3.0 source (from Ubuntu’s 2011 archive)
- `external/maloc-1.5/` — MALOC 1.5.0 source (from the official FETK distribution)

TAPBS applies four small patches to APBS 1.3 (in `src/patch_*`) to expose streamlined constructors, support precomputed maps, and reduce memory use.

The historical installer **`installTAPBS.sh`** is also included verbatim as `scripts/installTAPBS_2011.sh` for transparency. You don’t need to run it; the modern instructions below are safer.

---

## Relationship to Karlsberg 2.0

TAPBS computes the PB/LPBE electrostatics that Karlsberg uses to predict protonation and redox patterns in proteins.  
You can run TAPBS standalone or as part of the Karlsberg pipeline.

---

## Authors (original)

- **Gernot Kieseritzky** <gernotf@chemie.fu-berlin.de>  
- **Ilkay Sakalli** <sakalli@chemie.fu-berlin.de>  
Macromolecular Modelling Group (E.-W. Knapp), **Freie Universität Berlin**

---

## License

**GPL-2.0-or-later.** See `COPYING` for the full text.  
A short `LICENSE` stub is included for modern tooling; the canonical license remains `COPYING`.

Third-party licenses are preserved in their original distributions; see `THIRD_PARTY_NOTICES.md`.

---

## Documentation

Open `doc/userguide.html` in your browser for theory, parameters, and usage.  
Example runs are provided under `example/`.

---

## Building the historical stack (APBS 1.3 + MALOC 1.5)

> This reproduces the 2011 environment targeted by TAPBS 0.2.

### Prerequisites

**This script does not install OS packages for you.** Please pre-install dependencies:

- C/C++ toolchain (`gcc`, `g++`)
- Fortran compiler (`gfortran` recommended)
- BLAS (system `libblas` or `libopenblas`)
- `readline` development headers
- GNU Autotools (already generated `configure` is shipped)

You can use the commands:
- macOS (Homebrew):
  - `brew install gcc openblas`
- Linux (examples):
  - Debian/Ubuntu: `sudo apt-get install build-essential gfortran patch libopenblas-dev`
  - Fedora/RHEL: `sudo dnf install gcc gcc-c++ gcc-gfortran patch openblas-devel`
  - Arch: `sudo pacman -S base-devel gfortran openblas`

### Semi-Automated Build (modern archival setup)

Use the single cross-platform installer under `scripts/`:
- `scripts/install_tapbs_modern.sh` — builds MALOC 1.5, patches & builds APBS 1.3, then builds TAPBS 0.2. Installs into a local prefix without touching system paths.

Simply run:
- `bash scripts/install_tapbs_modern.sh`
- Optional custom prefix: `PREFIX=/your/path bash scripts/install_tapbs_modern.sh`


**Notes:**
> The script reproduces the exact autoconf cache hints and compiler flags required by the 2011 code, including disabling Python and tools in APBS and forcing OpenBLAS (`-lopenblas`).
> MALOC (the Minimal Abstraction Layer for Object-oriented C) is an APBS/FETk dependency, not the C allocator.  
  > The included historical version ensures compatibility with TAPBS’s APBS patches.

### Manual build (advanced, matches the script)

1) MALOC 1.5
- `cd external/maloc-1.5`
- `./configure --prefix="<prefix>"`
- `make -j && make install`

2) APBS 1.3 (patched)
- `cd external/apbs-1.3-source`
- `patch src/mg/vpmg.c            ../../src/patch_vpmgc`
- `patch src/mg/apbs/vpmg.h       ../../src/patch_vpmgh`
- `patch src/generic/vpbe.c       ../../src/patch_vpbec`
- `patch src/generic/apbs/vpbe.h  ../../src/patch_vpbeh`
- Environment (proven): set `CC`, `CXX`, `F77`; use `FFLAGS="-fallow-argument-mismatch -fno-second-underscore"`, and link with OpenBLAS.
- Configure (exactly):
  - `ac_cv_f77_compiler_gnu=yes ac_cv_f77_mangling="lower case, underscore, no extra underscore" ac_cv_lib_blas___dscal=yes ./configure --prefix="<prefix>" --with-blas='-lopenblas' --disable-openmp --disable-zlib --with-python=no --disable-tools`
- Build & install:
  - `make -j`
  - `make py_path=: install`   (forces Python parts to no-op)

3) TAPBS 0.2
- From repo root:
  - `CPPFLAGS="-include unistd.h -include sys/times.h" CFLAGS="-include unistd.h -include sys/times.h" CXXFLAGS="-include unistd.h -include sys/times.h" ./configure --prefix="<prefix>/tapbs" --with-apbs="<prefix>" --with-blas='-lopenblas -lgfortran'`
  - `make -j && make install`


### Cleaning and rebuilding

```bash
make clean
make distclean   # removes configuration
autoreconf -fi   # only if you edit configure.in or Makefile.am
```

- Then run the installer again or re-`./configure` with the same flags.

---

## Provenance

Archival source recovered from the Internet Archive (Wayback Machine).  
Original TAPBS 0.2 archive retrieved on **2025-10-08**.  

Third-party dependencies:
- **APBS 1.3.0** — obtained from [Ubuntu Launchpad source archive](https://launchpad.net/ubuntu/+source/apbs/1.3.0-2)
- **MALOC 1.5.0** — obtained from the [official FETK project site](http://www.fetk.org/codes/maloc.html)

For legal attributions and source origins, see `THIRD_PARTY_NOTICES.md`.

> **Transparency note:**  
> The historical installer script `installTAPBS.sh` is included **verbatim** for reference only.  
> It is *not* required for the modern, non-root build described above.

---

## Citing

Please cite the methodology referenced in the manual under section 6. Copyright.  
A machine-readable `CITATION.cff` is provided.

---

## Maintainer (archival)

**Gabriel Häusler**  
Department for Theoretical Biophysics, Institute for Theoretical Physics  
Johannes Kepler University Linz, Austria  
📧 gabriel.haeusler@posteo.com

Pull requests welcome for **build fixes** and **portability**.  
Algorithmic behavior should remain unchanged unless clearly documented.
