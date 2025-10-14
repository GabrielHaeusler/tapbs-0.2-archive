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

# Find exact compiler binaries and export them.
find_and_export_toolchain() {
  local os ver_list=("15" "14" "13" "12" "11" "10" "9")
  os="$(uname -s)"

  if [[ "$os" == "Darwin" ]]; then
    command -v brew >/dev/null 2>&1 || die "Homebrew is required. Install it first."
    local gcc_prefix
    gcc_prefix="$(brew --prefix gcc 2>/dev/null || true)"
    [[ -n "$gcc_prefix" && -d "$gcc_prefix/bin" ]] || die "Missing dependency: brew install gcc"

    local gc="" gp="" gf=""
    for v in "${ver_list[@]}"; do
      [[ -x "$gcc_prefix/bin/gcc-$v"      ]] && gc="$gcc_prefix/bin/gcc-$v"
      [[ -x "$gcc_prefix/bin/g++-$v"      ]] && gp="$gcc_prefix/bin/g++-$v"
      [[ -x "$gcc_prefix/bin/gfortran-$v" ]] && gf="$gcc_prefix/bin/gfortran-$v"
      if [[ -n "$gc" && -n "$gp" && -n "$gf" ]]; then break; fi
    done
    [[ -n "$gc" && -n "$gp" && -n "$gf" ]] || die "Homebrew gcc not found (looked for gcc-N, g++-N, gfortran-N). Try: brew install gcc"

    export CC="$gc"
    export CXX="$gp"
    export F77="$gf"
  else
    # Linux (and other UNIX): prefer versioned first, then unversioned.
    local gc="" gp="" gf=""
    for v in "${ver_list[@]}"; do
      command -v "gcc-$v"      >/dev/null 2>&1 && { gc="$(command -v gcc-$v)"; break; }
    done
    [[ -z "$gc" ]] && gc="$(command -v gcc || true)"
    for v in "${ver_list[@]}"; do
      command -v "g++-$v"      >/dev/null 2>&1 && { gp="$(command -v g++-$v)"; break; }
    done
    [[ -z "$gp" ]] && gp="$(command -v g++ || true)"
    for v in "${ver_list[@]}"; do
      command -v "gfortran-$v" >/dev/null 2>&1 && { gf="$(command -v gfortran-$v)"; break; }
    done
    [[ -z "$gf" ]] && gf="$(command -v gfortran || true)"

    [[ -n "$gc" && -n "$gp" && -n "$gf" ]] || die "Could not find GCC/g++/gfortran toolchain on PATH."
    export CC="$gc"
    export CXX="$gp"
    export F77="$gf"
  fi

  say "Using CC=$CC"
  say "Using CXX=$CXX"
  say "Using F77=$F77"
}

# Detect OpenBLAS location (macOS=Homebrew; Linux=system unless OPENBLAS_PREFIX set)
configure_blas_flags() {
  local os; os="$(uname -s)"
  if [[ "$os" == "Darwin" ]]; then
    local obp
    obp="$(brew --prefix openblas 2>/dev/null || true)"
    [[ -n "$obp" ]] || die "Missing dependency: brew install openblas"
    export LDFLAGS="-L${obp}/lib"
  else
    # On Linux, system link paths typically find -lopenblas without -L/-I.
    # Respect OPENBLAS_PREFIX if the user set it.
    if [[ -n "${OPENBLAS_PREFIX:-}" ]]; then
      export LDFLAGS="-L${OPENBLAS_PREFIX}/lib ${LDFLAGS:-}"
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
export FFLAGS="-fallow-argument-mismatch -fno-second-underscore"
export CFLAGS="-Wno-implicit-int -Wno-error=implicit-function-declaration"
# export CXXFLAGS="-O3"

# ===== MALOC 1.5 =====
say "==== MALOC 1.5 ===="
pushd external/maloc-1.5 >/dev/null
  ./configure --prefix="$PREFIX"
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
      --disable-tools

  make
  # Force python wrappers to no-op during install.
  make py_path=: install
popd >/dev/null

# ===== TAPBS 0.2 =====
say "==== TAPBS 0.2 ===="
TP_CPPFLAGS="-include unistd.h -include sys/times.h"
CPPFLAGS="$TP_CPPFLAGS" \
CFLAGS="$TP_CPPFLAGS" \
CXXFLAGS="$TP_CPPFLAGS" \
./configure \
    --prefix="$PREFIX/tapbs" \
    --with-apbs="$PREFIX" \
    --with-blas='-lopenblas -lgfortran'

make
make install

say "✅ Success. Binaries are under: $PREFIX/tapbs/bin"
if [[ "$OS" == "Darwin" ]]; then
  echo "Add to your shell:"
  echo "  export PATH=\"$PREFIX/tapbs/bin:\$PATH\""
  echo "  export DYLD_FALLBACK_LIBRARY_PATH=\"$PREFIX/lib:\${DYLD_FALLBACK_LIBRARY_PATH:-}\""
else
  echo "Add to your shell:"
  echo "  export PATH=\"$PREFIX/tapbs/bin:\$PATH\""
  echo "  export LD_LIBRARY_PATH=\"$PREFIX/lib:\${LD_LIBRARY_PATH:-}\""
fi
