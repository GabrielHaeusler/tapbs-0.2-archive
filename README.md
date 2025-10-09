# TAPBS 0.2

Front end to **APBS 1.3** that automates the many LPBE (linearized Poisson–Boltzmann) calculations needed for continuum-electrostatics **pKₐ** and **redox** workflows.
Prepares the inputs and aggregates the electrostatic energy terms consumed by **Karlsberg 2.0** (a.k.a. Karlsberg+).

> TAPBS supports local conformational changes, multiple charge states, regions with different dielectric constants, and implicit membrane slabs.
> See the original HTML manual in `doc/userguide.html`.

## Relationship to Karlsberg 2.0
TAPBS computes the PB/LPBE electrostatics that Karlsberg uses to predict protonation and redox patterns in proteins. You can run TAPBS standalone or as part of the Karlsberg pipeline.

## Authors
- **Gernot Kieseritzky** <gernotf@chemie.fu-berlin.de>
- **Ilkay Sakalli** <sakalli@chemie.fu-berlin.de>
Macromolecular Modelling Group (E.-W. Knapp), **Freie Universität Berlin**

## License
**GPL-2.0-or-later.** See `COPYING` for the full text.
A short `LICENSE` stub is included for modern tooling; the canonical license remains `COPYING`.

## Documentation
Open `doc/userguide.html` in your browser for theory, parameters, and usage.
Example runs are provided under `example/`.

## Requirements
- A C/C++/Fortran toolchain (GNU Autotools project).
- **APBS 1.3** available on your PATH (or configured via the TAPBS options described in the manual).

## Building (GNU Autotools)
```bash
./configure
make
# optional
sudo make install
```

To choose a different prefix:
```bash
./configure --prefix="$HOME/.local"
make && make install
```

If you are on Apple Silicon, you can build natively with Homebrew’s GCC:
```bash
brew install gcc
./configure CC=gcc-14
make
```

### Cleaning and rebuilding
```bash
make clean
make distclean # removes configuration
autoreconf -fi # (re)generate build system if you edit configure.in or Makefile.am
```

## Provenance
Archival source recovered from the Internet Archive (Wayback Machine).
Archived URL:
https://web.archive.org/web/20170628061521if_/http://agknapp.chemie.fu-berlin.de/karlsberg/tapbs-0.2.tar.gz
Retrieved: 2025-10-08

## Citing
Please cite the original TAPBS publications lsited in the manual under section 6. Copyright.
A machine-readable CITATION.cff is provided.

## Maintainers (archival)
**Gabriel Häusler**
Department for Theoretical Biophysics, Institute for Theoretical Physics
Johannes Kepler University Linz, Austria
📧 gabriel.haeusler@posteo.com
Pull requests welcome for build fixes and portability only; algorithmic behavior should remain unchanged unless clearly documented.