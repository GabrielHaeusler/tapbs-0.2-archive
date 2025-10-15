#!/usr/bin/env bash
# TAPBS 0.2 archival build (macOS + Linux)
# Builds MALOC 1.5, patches & builds APBS 1.3, then builds TAPBS 0.2.
# Installs to a local prefix; no root needed.
#
# Usage:
#   bash scripts/install_tapbs_archival.sh
# Optional:
#   PREFIX=/custom/path bash scripts/install_tapbs_archival.sh
#
set -euo pipefail

say() { printf "%b\n" "[$(date +%H:%M:%S)] $*"; }
die() { say "ERROR: $*"; exit 1; }

cores() {
  if command -v nproc >/dev/null 2>&1; then nproc
  elif [[ "$(uname -s)" == "Darwin" ]]; then sysctl -n hw.ncpu
  else echo 2
  fi
}

# Find exact compiler binaries and export them (portable).
find_and_export_toolchain() {
  local os; os="$(uname -s)"

  # Honor user overrides if already set (CC/CXX/FC/F77). Don't second-guess.
  if [[ -n "${CC:-}" || -n "${CXX:-}" || -n "${FC:-}" || -n "${F77:-}" ]]; then
    : "${CC:=cc}"
    : "${CXX:=c++}"
    : "${FC:=gfortran}"  # APBS 1.3 probes F77; set both.
    : "${F77:=gfortran}"
    say "Using user-specified toolchain: CC=$CC CXX=$CXX FC=$FC F77=$F77"
    return
  fi

  if [[ "$os" == "Darwin" ]]; then
    # macOS keeps Homebrew’s GCC versioned; keep your original Homebrew block.
    command -v brew >/dev/null 2>&1 || die "Homebrew is required. Install it first."
    local gcc_prefix gc gp gf
    gcc_prefix="$(brew --prefix gcc 2>/dev/null || true)"
    [[ -n "$gcc_prefix" && -d "$gcc_prefix/bin" ]] || die "Missing dependency: brew install gcc"
    for v in 15 14 13 12 11 10 9; do
      [[ -x "$gcc_prefix/bin/gcc-$v"      ]] && gc="$gcc_prefix/bin/gcc-$v"
      [[ -x "$gcc_prefix/bin/g++-$v"      ]] && gp="$gcc_prefix/bin/g++-$v"
      [[ -x "$gcc_prefix/bin/gfortran-$v" ]] && gf="$gcc_prefix/bin/gfortran-$v"
      [[ -n "$gc" && -n "$gp" && -n "$gf" ]] && break
    done
    [[ -n "$gc" && -n "$gp" && -n "$gf" ]] || die "brew gcc not found. Try: brew install gcc"
    export CC="$gc" CXX="$gp" FC="$gf" F77="$gf"
  else
    # Linux: prefer unversioned compilers per GNU/Autoconf conventions.
    # Fedora provides these via packages: gcc, gcc-c++, gcc-gfortran.
    export CC="${CC:-$(command -v gcc || command -v cc || true)}"
    export CXX="${CXX:-$(command -v g++ || command -v c++ || true)}"
    export FC="${FC:-$(command -v gfortran || true)}"
    export F77="${F77:-$FC}"
    [[ -n "$CC" && -n "$CXX" && -n "$F77" ]] || die "Missing GCC toolchain (install gcc gcc-c++ gcc-gfortran)."
  fi

  say "Using CC=$CC"
  say "Using CXX=$CXX"
  say "Using FC=$FC"
  say "Using F77=$F77"
}

# Detect OpenBLAS location (macOS=Homebrew; Linux=system unless OPENBLAS_PREFIX set)
configure_blas_flags() {
  local os; os="$(uname -s)"
  # Respect user-provided flags.
  if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists openblas; then
    local ob_cflags ob_libs
    ob_cflags="$(pkg-config --cflags openblas)"
    ob_libs="$(pkg-config --libs openblas)"
    export CPPFLAGS="${CPPFLAGS:-} ${ob_cflags}"
    export LDFLAGS="${LDFLAGS:-} ${ob_libs}"
    say "Found OpenBLAS via pkg-config."
  else
    if [[ "$os" == "Darwin" ]]; then
      local obp
      obp="$(brew --prefix openblas 2>/dev/null || true)" || true
      [[ -n "$obp" ]] || die "Missing dependency: brew install openblas"
      export LDFLAGS="-L${obp}/lib ${LDFLAGS:-}"
      export CPPFLAGS="-I${obp}/include ${CPPFLAGS:-}"
    else
      # Linux: usually -lopenblas works if openblas-devel is installed.
      # Respect OPENBLAS_PREFIX if user provided it.
      if [[ -n "${OPENBLAS_PREFIX:-}" ]]; then
        export LDFLAGS="-L${OPENBLAS_PREFIX}/lib ${LDFLAGS:-}"
        export CPPFLAGS="-I${OPENBLAS_PREFIX}/include ${CPPFLAGS:-}"
      fi
      # Don’t force anything else; configure step passes '-lopenblas'.
    fi
  fi
}

# ---------- main ----------
OS="$(uname -s)"
case "$OS" in
  Darwin)  say "Detected macOS.";  ;;
  Linux)   say "Detected Linux.";  ;;
  *)       die "Unsupported OS: $OS" ;;
esac

PREFIX="${PREFIX:-"$(pwd)/.local-prefix"}"
mkdir -p "$PREFIX"

find_and_export_toolchain
configure_blas_flags

# PARALLEL="$(cores)"
# say "Building with $PARALLEL parallel jobs."
say "Install prefix: $PREFIX"

# Historical flags that make the 2011 stack happy
export FFLAGS="-O3 -fallow-argument-mismatch -fno-second-underscore"
export CFLAGS="-O3 -Wno-implicit-int -Wno-error=implicit-function-declaration"
# export CXXFLAGS="-O3"

# ===== MALOC 1.5 =====
say "==== MALOC 1.5 ===="
pushd external/maloc-1.5 >/dev/null
  ./configure --prefix="$PREFIX" --enable-static --disable-shared
  make
  make install
popd >/dev/null

# ===== APBS 1.3 (patched) =====
say "==== APBS 1.3 (patched) ===="
pushd external/apbs-1.3-source >/dev/null
  # Apply TAPBS patches (paths match the original installer)
  patch src/mg/vpmg.c            ../../src/patch_vpmgc
  patch src/mg/apbs/vpmg.h       ../../src/patch_vpmgh
  patch src/generic/vpbe.c       ../../src/patch_vpbec
  patch src/generic/apbs/vpbe.h  ../../src/patch_vpbeh

  ac_cv_f77_compiler_gnu=yes \
  ac_cv_f77_mangling="lower case, underscore, no extra underscore" \
  ac_cv_lib_blas___dscal=yes \
  ./configure \
      --prefix="$PREFIX" \
      --with-blas='-lopenblas' \
      --disable-openmp \
      --disable-zlib \
      --with-python=no \
      --disable-tools \
      --enable-static \
      --disable-shared

  make
  # Force python wrappers to no-op during install.
  make py_path=: install
popd >/dev/null

# ===== TAPBS 0.2 =====

# Make configure respect --with-blas by replacing the hardcoded -lblas with ${BLAS}
# (configure sets BLAS from --with-blas earlier)
sed -i.bak '/APBS_LIBS="/ s/-lblas/${BLAS}/' configure

say "==== TAPBS 0.2 ===="
TP_CPPFLAGS="-include unistd.h -include sys/times.h"
CPPFLAGS="$TP_CPPFLAGS" \
CFLAGS="$TP_CPPFLAGS" \
CXXFLAGS="$TP_CPPFLAGS" \
./configure \
    --prefix="$PREFIX/tapbs" \
    --with-apbs="$PREFIX" \
    --with-blas='-lopenblas -lgfortran' \
    LIBS='-lopenblas -lgfortran'

make
make install

say "✅ Success. Binaries are under: $PREFIX/tapbs/bin"