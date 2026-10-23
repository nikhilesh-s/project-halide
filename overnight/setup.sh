#!/bin/bash
# One-time setup (needs internet, ~20-30 min, ~3 GB). Builds Quantum ESPRESSO 7.5 (pw.x only) from source
# inside overnight/.qe using a private conda toolchain, and downloads the exact pseudopotentials Anish used.
# There is no prebuilt QE for Apple Silicon on conda-forge or Homebrew. Touches nothing outside this folder.
set -euo pipefail
cd "$(dirname "$0")"
HERE="$PWD"
if [ ! -x .qe/bin/pw.x ]; then
  mkdir -p .mamba
  [ -x .mamba/bin/micromamba ] || curl -Ls https://micro.mamba.pm/api/micromamba/osx-arm64/latest | tar -xj -C .mamba bin/micromamba
  [ -d .qe/conda-meta ] || .mamba/bin/micromamba create -y -r "$HERE/.mamba/root" -p "$HERE/.qe" -c conda-forge \
      compilers openmpi fftw openblas make
  export PATH="$HERE/.qe/bin:$PATH"
  mkdir -p build && cd build
  [ -d q-e-qe-7.5 ] || curl -fsSL https://gitlab.com/QEF/q-e/-/archive/qe-7.5/q-e-qe-7.5.tar.gz | tar -xz
  cd q-e-qe-7.5
  ./configure MPIF90=mpif90 F90=gfortran CC=clang \
      BLAS_LIBS="-L$HERE/.qe/lib -lopenblas" LAPACK_LIBS="-L$HERE/.qe/lib -lopenblas" \
      FFT_LIBS="-L$HERE/.qe/lib -lfftw3" FFLAGS="-O2 -g0 -fallow-argument-mismatch" \
      LDFLAGS="-Wl,-rpath,$HERE/.qe/lib" > "$HERE/build/configure.log" 2>&1
  # configure does not detect conda's FFTW: select the FFTW3 backend and its headers by hand
  sed -i '' -e 's|^DFLAGS *=.*|& -D__FFTW3|' -e "s|^IFLAGS *=.*|& -I$HERE/.qe/include|" \
      -e "s|^FFT_LIBS *=.*|FFT_LIBS = -L$HERE/.qe/lib -lfftw3|" make.inc
  make -j"$(sysctl -n hw.ncpu)" pw > "$HERE/build/make.log" 2>&1
  cp bin/pw.x "$HERE/.qe/bin/pw.x"
  cd "$HERE"
fi
mkdir -p pseudo
while read -r md5 name; do
  [ -f "pseudo/$name" ] || curl -fsSL -o "pseudo/$name" "https://pseudopotentials.quantum-espresso.org/upf_files/$name"
  got=$(md5 -q "pseudo/$name")
  [ "$got" = "$md5" ] || { echo "MD5 MISMATCH for $name: $got (expected $md5)"; exit 1; }
done < pseudo_md5.txt
echo "SETUP OK: pw.x built, $(ls pseudo/*.UPF | wc -l | tr -d ' ') pseudopotentials verified"
