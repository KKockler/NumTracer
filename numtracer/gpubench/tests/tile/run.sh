#!/usr/bin/env bash
#
# run.sh -- regression tests for DiFfRG's divisibility-aware MDRange tile
#           (DIFFRG_DIVISIBLE_TILE, DiFfRG/common/kokkos.hh).
#
# Two tests, both needed. Either alone can pass for the wrong reason:
#
#   tiletest  unit-tests compute_divisible_tile() against tiles derived independently from the
#             Kokkos sources by ../../tools/tilewaste.py. Proves the chosen tile is right.
#
#   maptest   runs the real Integrator_p2_4D_3ang::map() over the production QCD_Nf2 shape
#             {64,32,6,6,6}, built twice (DIFFRG_DIVISIBLE_TILE=0 and =1), and requires the 64
#             results to be BIT-IDENTICAL -- re-tiling a parallel_for changes no arithmetic. It
#             also prints the tile each build actually launched with, which is the control that
#             stops "bit-identical" from silently meaning "the new path never fired".
#
# These are not in DiFfRG/tests/ because the DiFfRG test suite is not built here (header-only
# edits need only `cmake --install .`), and a test that never runs is not a test.
#
# Usage:  ./run.sh
# Env:    DIFFRG_SRC  DiFfRG include dir to test (default: the DiFfRG_KT working tree)
#         INSTALL     DiFfRG install prefix, for bundled Kokkos + link libraries
#         CCBIN       host compiler; CUDA 12.9 does not support the system GCC 16, so this
#                     defaults to g++-14
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DIFFRG_SRC=${DIFFRG_SRC:-/mnt/data/Documents/Uni/Code/DiFfRG_KT/DiFfRG/include}
INSTALL=${INSTALL:-$HOME/.local/share/DiFfRG}
CCBIN=${CCBIN:-/usr/bin/g++-14}
ARCH=${ARCH:-sm_89}
tmp=$(mktemp -d -t ntbtile_XXXX)
trap 'rm -rf "$tmp"' EXIT

wrapper="$INSTALL/bundled/bin/nvcc_wrapper"
common=(-I"$DIFFRG_SRC" -isystem "$INSTALL/bundled/include" -isystem "$INSTALL/include"
        -isystem /usr/include/suitesparse -isystem /opt/cuda/targets/x86_64-linux/include
        -std=c++20 -expt-relaxed-constexpr -extended-lambda -arch="$ARCH")

fail=0

# ---------------------------------------------------------------- tiletest (host-only, no link)
echo "== tiletest: compute_divisible_tile unit cases"
NVCC_WRAPPER_DEFAULT_COMPILER=$CCBIN "$wrapper" -ccbin "$CCBIN" "${common[@]}" -O2 \
  "$here/tiletest.cc" -o "$tmp/tiletest"
"$tmp/tiletest" || fail=1

# ---------------------------------------------------------------- maptest (needs the full link)
# The link line is long and install-specific; lift it from a built app rather than hardcoding it.
APP_LINK=${APP_LINK:-/mnt/data/Documents/Uni/Code/qcd-codes/vacuum/YangMills/SP/build/CMakeFiles/YangMills.dir/link.txt}
if [ ! -f "$APP_LINK" ]; then
  echo "== maptest: SKIPPED (no built app to lift link flags from; set APP_LINK)" >&2
  exit $fail
fi

linkflags=$(python3 - "$APP_LINK" <<'PY'
import shlex, sys, pathlib
toks = shlex.split(pathlib.Path(sys.argv[1]).read_text().strip())
out, drop = [], False
for t in toks:
    if drop: drop = False; continue
    if t == '-o': drop = True; continue
    if t.endswith('.o') or t.endswith('Flows.a'): continue
    if t.startswith('-Wl,--dependency-file'): continue
    if t.endswith('kokkos_launch_compiler') or t.endswith('nvcc_wrapper') or t.endswith('/c++'): continue
    out.append(t)
print(' '.join(shlex.quote(x) for x in out))
PY
)

echo
echo "== maptest: real map() over {64,32,6,6,6}, tile off vs on"
for v in 0 1; do
  NVCC_WRAPPER_DEFAULT_COMPILER=$CCBIN "$wrapper" -ccbin "$CCBIN" "${common[@]}" -O3 -DNDEBUG \
    -DDIFFRG_DIVISIBLE_TILE=$v "$here/maptest.cc" -o "$tmp/maptest_$v" $linkflags
done

# maptest reads parameter.json from its working directory
cd "$here"
for v in 0 1; do
  "$tmp/maptest_$v" > "$tmp/out_$v.txt" 2> "$tmp/tile_$v.txt"
  sed -n 's/^/   /p' "$tmp/tile_$v.txt"
done

# Control: the two builds must NOT have launched the same tile, or the comparison proves nothing.
if diff -q "$tmp/tile_0.txt" "$tmp/tile_1.txt" >/dev/null; then
  echo "  [FAIL] both builds launched the same tile - DIFFRG_DIVISIBLE_TILE had no effect"
  fail=1
elif ! grep -q 'waste=1.0000x' "$tmp/tile_1.txt"; then
  echo "  [FAIL] DIFFRG_DIVISIBLE_TILE=1 still launches masked lanes"
  fail=1
else
  echo "  [ok ] tile changed and the =1 build launches zero masked lanes"
fi

if diff -q "$tmp/out_0.txt" "$tmp/out_1.txt" >/dev/null; then
  echo "  [ok ] results bit-identical over all $(wc -l < "$tmp/out_0.txt") points"
else
  echo "  [FAIL] results differ:"
  diff "$tmp/out_0.txt" "$tmp/out_1.txt" | head -10
  fail=1
fi

# ------------------------------------------------- frameprobe: per-thread stack frame invariant
# The interpolators must reach the kernel BY REFERENCE. If someone reintroduces a by-value
# tuple_cat of (x, pos, m_args) in quadrature_integrator.hh, each thread gets a local-memory copy of
# every interpolator -- 3576 B on this 13-interpolator shape, which measured 2.3x of runtime on
# sm_89 (docs/NUMTRACER_PER_THREAD_FRAME.md). There is no compile switch for that any more, so this
# guards it as an absolute bound rather than an A/B: a correct build needs only a few tens of bytes
# of frame, a regressed one needs thousands.
FRAME_MAX_B=${FRAME_MAX_B:-256}
echo
echo "== frameprobe: per-thread stack frame (sm_80, sm_89)"
for arch in sm_80 sm_89; do
  NVCC_WRAPPER_DEFAULT_COMPILER=$CCBIN "$wrapper" -ccbin "$CCBIN" \
    -I"$DIFFRG_SRC" -isystem "$INSTALL/bundled/include" -isystem "$INSTALL/include" \
    -isystem /usr/include/suitesparse -isystem /opt/cuda/targets/x86_64-linux/include \
    -std=c++20 -expt-relaxed-constexpr -extended-lambda -arch=$arch -O3 -DNDEBUG \
    -Xptxas -v -c "$here/frameprobe.cc" -o "$tmp/frame_$arch.o" \
    2>&1 | grep -oE '[0-9]+ bytes stack frame' | grep -oE '^[0-9]+' | sort -rn > "$tmp/frames_$arch.txt"
  maxf=$(head -1 "$tmp/frames_$arch.txt")
  if [ "${maxf:-999999}" -le "$FRAME_MAX_B" ]; then
    echo "  [ok ] $arch largest device stack frame ${maxf} B (<= $FRAME_MAX_B)"
  else
    echo "  [FAIL] $arch largest device stack frame ${maxf} B exceeds $FRAME_MAX_B"
    echo "         a by-value interpolator copy is back in the map() path"
    fail=1
  fi
done

# The probe must still produce sane numbers, not just a small frame.
NVCC_WRAPPER_DEFAULT_COMPILER=$CCBIN "$wrapper" -ccbin "$CCBIN" "${common[@]}" -O3 -DNDEBUG \
  "$here/frameprobe.cc" -o "$tmp/frameprobe" $linkflags
"$tmp/frameprobe" > "$tmp/fr.txt"
if [ -s "$tmp/fr.txt" ] && ! grep -qiE 'nan|inf' "$tmp/fr.txt"; then
  echo "  [ok ] probe produced $(wc -l < "$tmp/fr.txt") finite results"
else
  echo "  [FAIL] probe output missing or non-finite"
  fail=1
fi

echo
[ $fail -eq 0 ] && echo "PASS" || echo "FAILURES"
exit $fail
