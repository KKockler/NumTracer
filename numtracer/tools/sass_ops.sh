#!/usr/bin/env bash
# sass_ops.sh — the GPU grading tool: fp64 instruction census of a CUDA TU on sm_89.
#
# Compiles <file.cu> to a cubin with the SAME production-relevant knobs the integrators use
# (NOTE __launch_bounds__: kernels without it measured ~14% slower — tests/gpu/README.md:78-81.
# The launch bounds live in the SOURCE (integrator.cuh pins 128); this script only reports
# whether any kernel in the TU carries them, so a bare probe kernel is flagged, not silently
# mis-graded) and reports:
#   - per-kernel registers / spill / stack from ptxas -v
#   - the fp64 opcode census from cuobjdump -sass: DADD DMUL DFMA DSETP MUFU.* and totals
#   - transcendental (MUFU) breakdown
#
# The fp64 instruction count is the success metric below the ~500-instr noinline cliff;
# registers/spill are NOT (that regime is settled — README:93-140). Above the cliff grade on
# wall-clock instead (the kernel goes memory-bound).
#
# Usage:  sass_ops.sh <file.cu> [extra nvcc flags...]
#         SASS_ARCH=sm_89 SASS_KEEP=1 sass_ops.sh probe.cu -I../tests/gpu
set -euo pipefail

[ $# -ge 1 ] || { echo "usage: $0 <file.cu> [nvcc flags...]" >&2; exit 1; }
src=$1; shift
arch=${SASS_ARCH:-sm_89}
out=$(mktemp -t sassops_XXXX.cubin)
trap '[ "${SASS_KEEP:-0}" = 1 ] || rm -f "$out"' EXIT

echo "=== nvcc -arch=$arch -O3 --expt-relaxed-constexpr -cubin $src $* ==="
nvcc -arch="$arch" -O3 --expt-relaxed-constexpr -cubin "$src" -o "$out" -Xptxas -v "$@" 2>&1 \
  | grep -E "Function properties|Used [0-9]+ registers|spill|bytes stack|ptxas info.*Compiling" \
  | sed 's/^ptxas info *: //'

echo
# launch bounds live in headers (integrator.cuh pins 128), so check the whole include closure
deps=$(nvcc -arch="$arch" --expt-relaxed-constexpr -M "$src" "$@" 2>/dev/null | tr ' \\' '\n' | grep -E '\.(cu|cuh|hh|hpp|h)$' | sort -u || true)
if grep -lq "__launch_bounds__\|LaunchBounds" $src $deps 2>/dev/null; then
  echo "launch bounds: present in include closure"
else
  echo "launch bounds: NOT FOUND — timing-relevant results need __launch_bounds__(128) (README:78-81)"
fi

echo
echo "=== fp64 opcode census (cuobjdump -sass) ==="
sass=$(cuobjdump -sass "$out")

# Per-kernel section split: cuobjdump prints 'Function : <mangled>' headers.
echo "$sass" | awk '
  /Function :/ { fn=$3; next }
  {
    # opcode is the first token after the /*addr*/ comment field
    for (i=1; i<=NF; i++) if ($i ~ /^(DADD|DMUL|DFMA|DSETP|MUFU|DADD32I|F2F|I2F.*F64|FSEL)/) {
      op=$i; sub(/;$/,"",op)
      split(op, parts, "."); base=parts[1]
      key = base (base=="MUFU" ? "."parts[2] : "")
      tot[key]++
      if (base ~ /^D/ || base=="MUFU") { fp64[fn]++; g64++ }
      break
    }
  }
  END {
    for (k in tot) printf "OP %-14s %d\n", k, tot[k]
    printf "TOTAL_fp64 %d\n", g64
    for (f in fp64) printf "KERNEL %s %d\n", f, fp64[f]
  }' > "$out.census"

grep '^OP '     "$out.census" | sort -k3 -n -r | awk '{printf "  %-14s %s\n", $2, $3}'
grep '^TOTAL'   "$out.census" | awk '{printf "  %-14s %s\n", "TOTAL fp64", $2}'
echo "  per-kernel fp64 totals:"
grep '^KERNEL ' "$out.census" | sort -k3 -n -r | awk '{printf "    %-70s %s\n", $2, $3}'
rm -f "$out.census"

[ "${SASS_KEEP:-0}" = 1 ] && echo "cubin kept at: $out"
exit 0
