#!/usr/bin/env bash
# sass.sh — L2/L2b: cross-architecture SASS census + occupancy of a flow's emitted kernels.
#
# Why this exists separately from ../../tools/sass_ops.sh: that script grades a hand-written probe
# on sm_89 and prints for humans. This one drives the whole production flow set across sm_80 (A100),
# sm_89 (RTX 4070) and sm_90 (H100) and emits CSV, because the whole point of this round is that the
# answer DIFFERS by architecture — the repo's "registers and occupancy are not the bottleneck"
# conclusion was established at 1:64 fp64 and does not transfer to the 1:2 datacenter parts.
#
# ptxas and cuobjdump cross-target happily, so every number here is obtainable without owning an
# A100 or an H100. Runtime (L3) still needs the real part.
#
# Usage:
#   sass.sh <flow_dir> [<flow_dir> ...]                 census, all three arches
#   NTB_ARCHES="sm_89" sass.sh <flow_dir>               one arch
#   NTB_EXTRA="-maxrregcount=128" sass.sh <flow_dir>    knob sweeps
#   NTB_CSV=results/l2.csv sass.sh <flow_dir>...        append machine-readable rows
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
arches=${NTB_ARCHES:-"sm_80 sm_89 sm_90"}
extra=${NTB_EXTRA:-}
csv=${NTB_CSV:-}
tmp=$(mktemp -d -t ntbsass_XXXX)
trap 'rm -rf "$tmp"' EXIT

# Registers per SM and max threads/SM, per architecture — for the occupancy model.
# (A100 sm_80 and H100 sm_90 both have 65536 regs/SM and 2048 max threads/SM; Ada sm_89 has
#  65536 regs/SM but only 1536 max threads/SM.)
regs_per_sm() { echo 65536; }
max_threads_per_sm() { case "$1" in sm_89) echo 1536;; *) echo 2048;; esac; }

if [ -n "$csv" ] && [ ! -s "$csv" ]; then
  mkdir -p "$(dirname "$csv")"
  echo "flow,arch,extra,regs,stack,spill_store,spill_load,dadd,dmul,dfma,dsetp,mufu,fp64_total,occ_pct,waves" > "$csv"
fi

printf "%-12s %-7s %6s %6s %7s %7s %8s %8s %8s %7s %7s %10s %7s\n" \
  flow arch regs stack spillS spillL DADD DMUL DFMA DSETP MUFU fp64_tot occ%

for flow in "$@"; do
  [ -f "$flow/kernels.hh" ] || continue
  name=$(basename "$flow")
  probe="$tmp/$name.cu"
  python3 "$here/mkprobe.py" "$flow" -o "$probe" >/dev/null || continue

  for arch in $arches; do
    cubin="$tmp/$name.$arch.cubin"
    ptxlog="$tmp/$name.$arch.ptxas"
    if ! nvcc -arch="$arch" -O3 --expt-relaxed-constexpr -cubin "$probe" -o "$cubin" \
              -Xptxas -v $extra > "$ptxlog" 2>&1; then
      echo "  $name $arch: COMPILE FAILED (see $ptxlog)" >&2
      sed -n '1,12p' "$ptxlog" >&2
      continue
    fi

    regs=$(grep -oE 'Used [0-9]+ registers' "$ptxlog" | grep -oE '[0-9]+' | head -1); regs=${regs:-0}
    stack=$(grep -oE '[0-9]+ bytes stack frame' "$ptxlog" | grep -oE '^[0-9]+' | head -1); stack=${stack:-0}
    sst=$(grep -oE '[0-9]+ bytes spill stores' "$ptxlog" | grep -oE '^[0-9]+' | head -1); sst=${sst:-0}
    sld=$(grep -oE '[0-9]+ bytes spill loads' "$ptxlog" | grep -oE '^[0-9]+' | head -1); sld=${sld:-0}

    # opcode census: first token of the instruction, after the /*addr*/ field
    read -r dadd dmul dfma dsetp mufu <<EOF
$(cuobjdump -sass "$cubin" | awk '
  { for (i=1;i<=NF;i++) if ($i ~ /^(DADD|DMUL|DFMA|DSETP|MUFU)(\.|;|$)/) {
      op=$i; sub(/;$/,"",op); split(op,a,"."); c[a[1]]++; break } }
  END { printf "%d %d %d %d %d", c["DADD"], c["DMUL"], c["DFMA"], c["DSETP"], c["MUFU"] }')
EOF
    fp64=$(( dadd + dmul + dfma + dsetp + mufu ))

    # theoretical occupancy at the production block size of 128 threads
    rsm=$(regs_per_sm); mts=$(max_threads_per_sm "$arch")
    if [ "$regs" -gt 0 ]; then
      # registers are allocated per warp in units of 256 on all these arches
      per_warp=$(( ((regs * 32 + 255) / 256) * 256 ))
      warps=$(( rsm / per_warp ))
      max_warps=$(( mts / 32 ))
      [ "$warps" -gt "$max_warps" ] && warps=$max_warps
      occ=$(( warps * 100 / max_warps ))
    else
      warps=0; occ=0
    fi

    printf "%-12s %-7s %6s %6s %7s %7s %8s %8s %8s %7s %7s %10s %6s%%\n" \
      "$name" "$arch" "$regs" "$stack" "$sst" "$sld" "$dadd" "$dmul" "$dfma" "$dsetp" "$mufu" "$fp64" "$occ"

    [ -n "$csv" ] && echo "$name,$arch,\"$extra\",$regs,$stack,$sst,$sld,$dadd,$dmul,$dfma,$dsetp,$mufu,$fp64,$occ,$warps" >> "$csv"
  done
done

exit 0
