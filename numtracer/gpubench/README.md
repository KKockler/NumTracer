# gpubench — GPU measurement harness for production NumTracer kernels

Measures the kernels the two production fRG codes actually run, without touching either code.
Everything is driven off **vendored copies** under `kernels/`, so an experiment can never disturb
the tree it came from. Results and conclusions live in [FINDINGS.md](FINDINGS.md).

Target hardware is **A100 (sm_80) and H100 (sm_90)**, with the local RTX 4070 Laptop (sm_89) as a
secondary target. That distinction matters: fp64 is 1:64 on Ada and 1:2 on the datacenter parts, so
conclusions do not transfer between them and every layer is parameterised by architecture.

## Layers

| layer | tool | measures | needs the target GPU? |
|---|---|---|---|
| L1 static | `tools/census.py` | fp64 SSA ops, transcendental total vs distinct, interpolator count, per-function size histogram | no |
| L2 SASS | `tools/sass.sh` | regs / stack / spill, DADD·DMUL·DFMA·DSETP·MUFU, for sm_80 **and** sm_89 **and** sm_90 | no |
| L2b occupancy | folded into `sass.sh` | theoretical occupancy & warps/SM at the production block size | no |
| L3 runtime | `harness/bench.cc` | ms per `map()`, median-of-N + spread, plus a checksum | yes |
| L4 profile | `nsys` / `ncu` | phase-1 vs phase-2 split, stall reasons | yes |

`ptxas` and `cuobjdump` cross-target happily, so **L1/L2/L2b give real A100/H100 answers from a
machine that has neither**. Only L3/L4 need the actual part.

**Grading rule.** On sm_89 the fp64 instruction count is the metric. On sm_80/sm_90 it is not
sufficient on its own — occupancy and spill matter there — so L2b and L4 decide. On the laptop, a
runtime delta smaller than the reported min/max spread is not a result; use L2.

## Quick start

```bash
./vendor.sh                                     # refresh copies from the QCD codes

# L1 — no build required
python3 tools/census.py --code nf2 kernels/nf2/flows/*/

# L2/L2b — cross-architecture, no GPU required
NTB_CSV=results/l2.csv ./tools/sass.sh kernels/nf2/flows/*/
NTB_ARCHES=sm_90 NTB_EXTRA=-maxrregcount=128 ./tools/sass.sh kernels/nf2/flows/ZA4

# L3 — runtime on the local GPU
cmake -S . -B build/baseline -DNTB_FLOWS=kernels/ym -DCMAKE_BUILD_TYPE=Release
cmake --build build/baseline -j8
./build/baseline/ntb_bench -p kernels/ym/parameter.toml
```

## Variant rewrites

Each tool rewrites a vendored kernel into `variants/<name>/` so a lever can be sized in minutes,
before anything in the emitter or in DiFfRG is changed. All of them preserve semantics and are
gated on reproducing the baseline checksum bit-for-bit.

| tool | lever | verdict |
|---|---|---|
| `tools/interpshare.py` | hoist one interpolator `index()` per distinct argument | **1.13–1.26x measured** |
| `tools/noinline.py` | apply the device `__noinline__` size gate at a chosen threshold | open; threshold is arch-dependent |
| `tools/fillcse.py` | CSE `fill()` | inert — nvcc already does it |
| `tools/mkprobe.py` | build a standalone CUDA probe of the NumTracer emission alone | used by `sass.sh` |

```bash
cp -a kernels/ym variants/interpshare
for f in variants/interpshare/flows/*/; do
  python3 tools/interpshare.py kernels/ym/flows/$(basename $f)/kernel.hh -o $f/kernel.hh
done
cmake -S . -B build/interpshare -DNTB_FLOWS=variants/interpshare -DCMAKE_BUILD_TYPE=Release
cmake --build build/interpshare -j8 && ./build/interpshare/ntb_bench -p kernels/ym/parameter.toml
```

The bench prints a checksum per flow. **A semantics-preserving variant must reproduce the baseline
checksum exactly**; any deviation is a bug in the rewrite, not a tradeoff.

## Scope limits

- Both flow sets are built and timed. `CMakeLists.txt` picks the driver from the vendored tree
  (`harness/bench_nf2.cc` when `flows/Zq` exists, else `harness/bench.cc`).
- **Timing discipline.** The GPU idles at 210 MHz and needs ~1 s of load to reach boost, so both
  drivers burn a fixed 2.5 s (`NTB_BURN_MS`) before timing. Alternate variant order, run >= 2 rounds,
  discard round 1, compare ratios not absolute ms — and prefer the fp64 SASS count for the long Nf2
  flows, where sustained load throttles hard enough to invert a wall-clock comparison.
- Interpolator fixtures use a smooth analytic profile, not physical dressings. Correct for timing
  and for A/B comparison; it is **not** physics and nothing here claims otherwise.
- `numtracer_add_numtrace()` is deliberately not called: the harness measures committed kernels and
  must never trigger a regeneration as a build side effect.
