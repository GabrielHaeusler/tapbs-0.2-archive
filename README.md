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

The historical installer **`installTAPBS.sh`** is also included verbatim for transparency. You don’t need to run it; the modern instructions below are safer.

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

- C/C++ toolchain (`gcc`, `g++`)
- Fortran compiler (`gfortran` recommended)
- BLAS (system `libblas` or `libopenblas`)
- `readline` development headers
- GNU Autotools (already generated `configure` is shipped)

### Install prefix

```bash
export PREFIX="$HOME/.local/tapbs-0.2"
mkdir -p "$PREFIX"
```

### 1) Build MALOC 1.5

```bash
cd external/maloc-1.5
./configure --prefix="$PREFIX"
make -j && make install
cd ../..
```

### 2) Patch and build APBS 1.3

```bash
cd external/apbs-1.3-source
patch -p0 < ../../src/patch_vpbec
patch -p0 < ../../src/patch_vpbeh
patch -p0 < ../../src/patch_vpmgc
patch -p0 < ../../src/patch_vpmgh
./configure --prefix="$PREFIX" --with-blas='-lblas' --disable-openmp --disable-zlib
make -j && make install
cd ../..
```

### 3) Build TAPBS 0.2

```bash
./configure --with-apbs="$PREFIX" --with-blas='-lblas -lgfortran'
make -j
make install   # optional; installs under $PREFIX if you add --prefix to configure
```

### Environment

```bash
echo 'export PATH="'$PREFIX'/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="'$PREFIX'/lib:${LD_LIBRARY_PATH:-}"' >> ~/.bashrc
# reload your shell or source ~/.bashrc
```

> **Note:** MALOC (the Minimal Abstraction Layer for Object-oriented C) is an APBS/FETk dependency, not the C allocator.  
> The included historical version ensures compatibility with TAPBS’s APBS patches.

---

## Cleaning and rebuilding

```bash
make clean
make distclean   # removes configuration
autoreconf -fi   # only if you edit configure.in or Makefile.am
```

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
