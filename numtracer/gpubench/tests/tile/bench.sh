#!/usr/bin/env bash
#
# bench.sh -- measure DIFFRG_DIVISIBLE_TILE on the LOCAL GPU.
#
# Alternates the two binaries round-robin.
#
# (The by-reference argument passing this used to A/B alongside the tile is no longer a switch: it
#  measured 2.323x on sm_89 and is now unconditional. See docs/NUMTRACER_PER_THREAD_FRAME.md.)
#
# Absolute times drift as the part heats, so only same-round ratios are meaningful -- and round 1 is discarded because the
# 4070 Laptop idles at 210 MHz (bench.cc burns for a fixed duration first, but the first process of
# a session still sees a colder card). This is the trap that produced a bogus "1.3-1.5x" tile
# result in docs/NUMTRACER_GPU_INVESTIGATION.md 5; do not shortcut it.
#
# Usage: ./bench.sh [rounds] [iters] [burn_seconds]
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
rounds=${1:-4}
iters=${2:-200}
burn=${3:-3.0}
DIFFRG_SRC=${DIFFRG_SRC:-/mnt/data/Documents/Uni/Code/DiFfRG_KT/DiFfRG/include}
INSTALL=${INSTALL:-$HOME/.local/share/DiFfRG}
CCBIN=${CCBIN:-/usr/bin/g++-14}
ARCH=${ARCH:-sm_89}
APP_LINK=${APP_LINK:-/mnt/data/Documents/Uni/Code/qcd-codes/vacuum/YangMills/SP/build/CMakeFiles/YangMills.dir/link.txt}
tmp=$(mktemp -d -t ntbbench_XXXX)
trap 'rm -rf "$tmp"' EXIT

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

echo "building tile off/on for $ARCH ..."
for tile in 0 1; do
  NVCC_WRAPPER_DEFAULT_COMPILER=$CCBIN "$INSTALL/bundled/bin/nvcc_wrapper" -ccbin "$CCBIN" \
    -I"$DIFFRG_SRC" -isystem "$INSTALL/bundled/include" -isystem "$INSTALL/include" \
    -isystem /usr/include/suitesparse -isystem /opt/cuda/targets/x86_64-linux/include \
    -std=c++20 -expt-relaxed-constexpr -extended-lambda -arch="$ARCH" -O3 -DNDEBUG \
    -DDIFFRG_DIVISIBLE_TILE=$tile \
    "$here/bench.cc" -o "$tmp/bench_t${tile}" $linkflags 2>&1 | grep -viE "warning|remark|^$" || true
done

cd "$here"   # bench.cc reads parameter.json from cwd
declare -A sum
echo
printf "%-7s %-16s %-16s\n" round "t0(kokkos tile)" "t1(divisible)"
for r in $(seq 1 "$rounds"); do
  line=$(printf "%-7s " "$r")
  for key in t0 t1; do
    out=$("$tmp/bench_$key" "$iters" "$burn")
    ms=$(echo "$out" | grep -oE 'per_iter_ms=[0-9.]+' | cut -d= -f2)
    cs=$(echo "$out" | grep -oE 'checksum=[-0-9.e+]+' | cut -d= -f2)
    echo "$r $key $ms $cs" >> "$tmp/raw.txt"
    line+=$(printf "%-16s " "$ms")
  done
  echo "$line"
done

echo
python3 - "$tmp/raw.txt" "$rounds" <<'PY'
import sys, statistics, collections
rows = [l.split() for l in open(sys.argv[1])]
rounds = int(sys.argv[2])
keep = {r for r in range(2, rounds + 1)} if rounds > 1 else {1}
by = collections.defaultdict(list)
cs = collections.defaultdict(set)
for rnd, key, ms, chk in rows:
    cs[key].add(chk)
    if int(rnd) in keep: by[key].append(float(ms))

label = {'t0': 'kokkos tile (control)', 't1': 'divisibility tile'}
if len(keep) == 1:
    print(f"NOTE: only round(s) {sorted(keep)} kept - run more rounds for a trustworthy median")
else:
    print(f"medians over rounds {sorted(keep)} (round 1 discarded as warm-up)")
base = statistics.median(by['t0'])
print(f"{'variant':<22s} {'ms/iter':>9s} {'spread':>8s} {'vs baseline':>12s}")
for k in ('t0', 't1'):
    v = by[k]
    m = statistics.median(v)
    print(f"{label[k]:<22s} {m:>9.4f} {max(v)-min(v):>8.4f} {base/m:>11.3f}x")

# re-tiling moves no arithmetic, so this must match exactly
allcs = {k: sorted(v) for k, v in cs.items()}
if allcs['t0'] == allcs['t1']:
    print("\nchecksum: identical (expected - re-tiling moves no arithmetic)")
else:
    print(f"\nWARNING: checksum differs between tilings: {allcs['t0']} vs {allcs['t1']}")
PY
