# Contributing (Archival)

Thanks for helping preserve TAPBS.

## Scope
This repository is an archival publication. Changes should be minimal and aimed at:
- Portability fixes (compiler warnings, 64-bit/ARM, Autotools refresh)
- Reproducible builds (no algorithmic changes by default)
- Documentation clarifications that do not alter scientific behavior

## Build system
TAPBS uses GNU Autotools. Typical workflow:

./configure
make
make install # optional

If `configure` fails to detect your host, refresh `config.guess`/`config.sub` or run `autoreconf -fi`.

## Testing
Include a small reproducer (an `example/` case and expected log/output).
Numerical results should remain stable within the tolerances discussed in the manual.

## License and authorship
- Original code is **GPL-2.0-or-later** (see `COPYING`).
- Preserve all copyright and attribution notices.
- Do not relicense or remove author headers.
