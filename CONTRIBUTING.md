# Contributing (Archival)

Thanks for helping preserve TAPBS 0.2.

This repository is an **archival re-release**, not an active development project.  
Changes should focus on:

- Build or portability fixes (compiler warnings, 64-bit/ARM)  
- Autotools maintenance (configure updates)  
- Documentation or reproducibility improvements  

Algorithmic or scientific behavior must not change without full documentation.

---

## Building locally

TAPBS uses **GNU Autotools** and expects its historical dependencies (APBS 1.3 + MALOC 1.5) located under `external/`.

```bash
./configure --with-apbs="$HOME/.local/tapbs-0.2" --with-blas='-lblas -lgfortran'
make -j
make install  # optional
```

If `configure` fails to detect your platform:
```bash
autoreconf -fi
```

---

## Testing

Use one of the example inputs under `example/`.  
Expected numerical values should match the reference outputs within floating-point tolerance.

---

## Licensing and attribution

- TAPBS 0.2 is **GPL-2.0-or-later** (see `COPYING`).  
- Preserve all copyright and authorship notices.  
- Do not relicense or remove author headers.  

Third-party code in `external/` retains its original licenses (see `THIRD_PARTY_NOTICES.md`).
