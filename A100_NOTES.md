# A100 measurement session — results log

Session run per `NUMTRACER_A100_SESSION.md`. **Measurement only** — no commits, no default changes.

## 0. Environment

| item | value |
|---|---|
| date | 2026-08-10 |
| cluster | Noctua2 (PC2 Paderborn) |
| node | `n2gpu1214` |
| SLURM job | 33988800 |
| GPU | NVIDIA A100-SXM4-40GB (sm_80), 40960 MiB |
| driver | 580.173.02 |
| ECC | Enabled |
| SM clock | max 1410 MHz (idle 210 MHz) |
| CPU | 64 cores allocated, 503 GB RAM |
| CUDA | 12.8.0 (nvcc V12.8.61) |
| GCC | 14.2.0 |
| CMake | 3.31.3 |
| Python (tools) | 3.13.1 — system `python3` is 3.6.8 and **cannot run `noinline.py`** (`from __future__ import annotations`); use `/opt/software/pc2/EB-SW/software/Python/3.13.1-GCCcore-14.2.0/bin/python3` |

Module set (`qcd-codes/dev_GPU.sh`, identical to `DiFfRG/dev_GPU.sh`):
```
devel/CMake/3.31.3-GCCcore-14.2.0
numlib/GSL/2.8-GCC-14.2.0
numlib/OpenBLAS/0.3.29-GCC-14.2.0
system/CUDA/12.8.0
```

### Paths (differ from the session doc's `$T` layout)

| tree | path |
|---|---|
| DiFfRG source | `~/scratch/Franz/DiFfRG` (branch `feature/fast_interpolators`, HEAD `7ab9cc5fa`) |
| DiFfRG build dir | `~/scratch/Franz/DiFfRG/bin` (library-only build; **not** `build_GPU`, which is the superbuild) |
| DiFfRG install | `~/.local/share/DiFfRG_GPU` |
| apps | `~/Code/qcd-codes/vacuum/{QCD_Nf2/no_mesons, YangMills/SP}` |
| tools + notes | `~/Code/numtracer/{gpubench-tools, ref-results}` |

DiFfRG install was already current at session start (rebuilt + installed 2026-08-10 10:27/10:30,
contains `fastInterpLookups`, `DIFFRG_LAUNCH_BOUNDS` default 128). Installed headers were
`touch`ed before the app builds (§11 trap) — both apps then recompiled from scratch.

## 1. Pre-flight findings

### 1.1 nf2 ZA4 numtrace regeneration is nondeterministic in `tr0`

`flows/ZA4/numtrace.json` ships with `"generated": 0`, so the first configure scheduled a
re-emit — and, unlike the session doc's expectation, it **succeeded**: the numtrace step
compiles and runs a C++ generator (`gen/gen_za4_qcd_num*.cpp`), no Wolfram involved.

The re-emitted `kernels.hh` differs from the committed one:

| | committed | re-emitted |
|---|---|---|
| functions | 75 | 75 |
| distinct numeric literals | 2006 | 2006 |
| `fma()` calls | 5692 | 5688 |
| functions with differing op count | — | 1 (`tr0`: 1312 → 1326 SSA) |
| diff line range | — | 31–1343; `tr0` spans lines 30–1357 |

`tr1`…`tr74` are byte-identical. Only `tr0` moves, and `tr0` is the only function large enough
to be split across the generator's parallel CSE workers (`tracing ZA4 (W=8, WB=48)`), whose
scheduling is not order-deterministic. Effect: +14 SSA ops = +1.1% on `tr0`, +0.13% on the flow.
Same math, same literal pool.

**Action taken:** restored the committed `kernels.hh` (`git checkout`, md5 `28d180cb…` verified)
and set `"generated": 1` in the manifest, so no variant tree can re-roll this. `numtrace.json`
is the only modified tracked file in the working tree.

### 1.2 `__noinline__` size-gate coverage (static, from `noinline.py --report`)

Currently **0** `noinline` decorators in all 16 production flows, both codes.

QCD_Nf2 — flows that gate at all (all others 0 at every threshold):

| flow | total SSA | min=200 | min=300 | min=500 | min=800 |
|---|---|---|---|---|---|
| ZA4 | 10415 | 8 fns / 34% ops | 7 fns / 31% | 1 fn / 12% | 1 fn / 12% |
| ZAAqbq2 | 3296 | 3 fns / 19% | 0 | 0 | 0 |
| ZAqbq4 | 4985 | 1 fn / 12% | 1 fn / 12% | 1 fn / 12% | 0 |

YangMills — only ZA4 gates:

| flow | total SSA | min=200 | min=300 | min=500 | min=800 |
|---|---|---|---|---|---|
| ZA4 | 2000 | 3 fns / 96% ops | 2 fns / 82% | 1 fn / 65% | 1 fn / 65% |

Note **800 ≠ off** for either code: `tr0` (nf2 1312 SSA, YM 1302 SSA) still gates at 800.
A true "off" arm is the unmodified tree.

---

## 2. Experiment 1 — baseline + `fastInterpLookups` A/B

Binaries: `no_mesons/build/QCD`, `SP/build/YangMills`, both stock (LB 128, no noinline),
`cuobjdump -lelf` confirms **sm_80** only. Configs: nf2 `run_gate/parameter.json` (OFF) vs
`run_fastinterp/parameter.json` (ON); YM `run_gate/parameter.toml` (OFF) vs
`run_fastinterp/parameter.toml` (ON), the two differing in exactly one line
(`fastInterpLookups`) — YM's committed `parameter.toml` ships with it **true**.

### 2.1 Physics gate — both codes, both flag states: PASS

| code | flag | m2A | target | rel. dev | Zc(0) | rel. dev |
|---|---|---|---|---|---|---|
| Nf2 | OFF | 1.1731490534726379 | 1.1731490534726277 | 8.7e-15 | 0.6888524203451547 | 3.0e-15 |
| Nf2 | ON  | 1.1731490534727145 | " | 7.4e-14 | 0.6888524203451604 | 1.1e-14 |
| YM  | OFF | 6037.244194877355 | 6037.244194877355 | **0** | 0.9316318995037499 | 6.4e-16 |
| YM  | ON  | 6037.244194877391 | " | 6.0e-15 | 0.9316318995037505 | **0** |

All ≤1e-8 by twelve orders of magnitude. YM OFF reproduces the reference m2A **bit-identically**.
Workloads are fixed as designed: Nf2 350 steps every run, YM 240 steps / 425 RHS every run.

### 2.2 Timing — 4 alternated rounds, round 1 discarded as warmup

Nf2 (`wall` seconds):

| round | OFF | ON | OFF/ON |
|---|---|---|---|
| 1 (warmup) | 22.72 | 23.62 | 0.9619 |
| 2 | 22.34 | 23.96 | 0.9324 |
| 3 | 22.67 | 23.72 | 0.9557 |
| 4 | 22.52 | 23.63 | 0.9530 |
| **median (2–4)** | | | **0.9530** |

YM (`wall` seconds):

| round | OFF | ON | OFF/ON |
|---|---|---|---|
| 1 (warmup) | 5.29 | 5.16 | 1.0252 |
| 2 | 5.10 | 5.09 | 1.0020 |
| 3 | 5.08 | 5.15 | 0.9864 |
| 4 | 5.07 | 5.14 | 0.9864 |
| **median (2–4)** | | | **0.9864** |

### 2.3 Result — the sign flips relative to the 4070

| code | 4070 (ON gain) | **A100 (ON gain)** |
|---|---|---|
| YangMills | 1.13x | **0.986x** (ON costs 1.4%) |
| QCD_Nf2 | 1.02x | **0.953x** (ON costs 4.9%) |

`fastInterpLookups` is **neutral-to-harmful on A100**: the tabulated inverse coordinate transform
trades `log`/`log1p`/`asinh` arithmetic for a table load. On the 1:64-fp64 Ada part that was a
clear win; on the 1:2 A100 the transcendental is cheap and the extra memory traffic — into kernels
already at 12.5% occupancy, i.e. with almost no warps to hide latency — is a net loss. Nf2 loses
more than YM, consistent with Nf2 being the more interpolator-heavy and more occupancy-starved code.

**Deliverable:** the OFF state is the A100 baseline for everything below. Note YM's committed
`parameter.toml` ships `fastInterpLookups = true`, which on A100 costs ~1.4%.

### 2.4 A100 baseline walls (denominator for the rest of the session)

| code | A100 wall (flag OFF) | 4070 reference | speedup |
|---|---|---|---|
| QCD_Nf2 (`run_gate`, 350 steps) | **22.5 s** | 106 s | 4.7x |
| YangMills (`parameter.toml`, 240 steps / 425 RHS) | **5.08 s** | 22 s | 4.3x |

## 3. Experiment 2 — per-flow GPU profile

`nsys profile` on the stock binaries. Three profiles: nf2 OFF, nf2 ON, YM OFF.
(The doc suggested profiling the flag-ON state as "the intended production state"; §2.3 shows
OFF is the better state on A100, so OFF is what is profiled as the baseline — ON is profiled too,
to locate the 4.9%.)

Instance counts confirm the workloads: nf2 **714 RHS evaluations** per solve, exactly the 4070
figure; YM 989/2 ≈ 494 launches of the big flows over 425 RHS.

### 3.1 QCD_Nf2, `fastInterpLookups` OFF — total GPU kernel time 14.512 s

| flow | GPU time (s) | share (A100) | share (4070) | mean per launch (µs) |
|---|---|---|---|---|
| ZA4 | 4.748 | **32.7%** | 39.1% | 3322.9 |
| ZAAqbq2 | 4.040 | **27.8%** | 26.9% | 2827.4 |
| ZAAqbq1 | 3.182 | **21.9%** | 17.1% | 2226.7 |
| ZAqbq4 | 0.422 | 2.9% | — | 295.6 |
| ZAqbq7 | 0.389 | 2.7% | — | 272.2 |
| ZAqbq1 | 0.382 | 2.6% | — | 267.4 |
| ZA | 0.313 | 2.2% | — | 34.6 |
| Zc | 0.293 | 2.0% | — | 32.4 |
| Zq | 0.261 | 1.8% | — | 34.6 |
| ZA3 | 0.257 | 1.8% | — | 179.5 |
| ZAcbc | 0.223 | 1.5% | — | 156.2 |

Top three = **82.4%** (4070: 83%) — but the mix inside those three has shifted: ZA4, the
77%-NumTracer-arithmetic kernel, drops 39.1% → 32.7%, while the 57–58%-interpolator pair
ZAAqbq1+2 rises 44.0% → 49.7%. That is exactly the expected consequence of fp64 going 1:64 → 1:2:
the polynomial-trace arithmetic gets ~40x cheaper relative to memory, so the interpolator-bound
kernels become the larger share. The ZAqbq family holds at ~8% (4070: 9.3%).

### 3.2 QCD_Nf2, `fastInterpLookups` ON — total GPU kernel time 15.676 s (+8.0%)

| flow | OFF (s) | ON (s) | Δ |
|---|---|---|---|
| ZA4 | 4.748 | 5.010 | +5.5% |
| ZAAqbq2 | 4.040 | 4.455 | +10.3% |
| ZAAqbq1 | 3.182 | 3.499 | +10.0% |
| ZAqbq4 | 0.422 | 0.452 | +7.1% |
| ZA | 0.313 | 0.336 | +7.3% |
| **total** | **14.512** | **15.676** | **+8.0%** |

The table transform costs GPU time in *every* flow, worst in the two interpolator-dominated ones
(+10%) — i.e. the tabulated inverse transform is slower than the `log`/`log1p`/`asinh` it replaces
on this part, not merely a wash. This is the mechanism behind the §2.3 wall-clock sign flip.

### 3.3 YangMills, `fastInterpLookups` OFF — total GPU kernel time 2.434 s

| flow | GPU time (s) | share (A100) | share (4070) | mean per launch (µs) |
|---|---|---|---|---|
| ZA4 | 2.017 | **82.9%** | 86.9% | 2039.2 |
| ZAcbc | 0.163 | 6.7% | — | 164.6 |
| ZA3 | 0.157 | 6.4% | — | 158.3 |
| ZA | 0.050 | 2.1% | — | 27.1 |
| Zc | 0.048 | 2.0% | — | 25.9 |

### 3.4 GPU time vs wall — a large and growing CPU share

Against the unprofiled baseline walls of §2.4:

| code | GPU kernel time | wall | **GPU share** | non-GPU |
|---|---|---|---|---|
| QCD_Nf2 | 14.512 s | 22.5 s | **64.5%** | 8.0 s (35.5%) |
| YangMills | 2.434 s | 5.08 s | **47.9%** | 2.6 s (52.1%) |

This is the most consequential number of the session and it is not in the 4070 picture at all.
On A100 the kernels got ~4.5x faster while the host-side work did not, so **over half of the
YangMills solve and over a third of the Nf2 solve is now outside the GPU kernels** (host
timestepper, launch latency, per-step CPU work). Every GPU-side optimisation below is bounded by
those shares — a 10% kernel win is 6.5% e2e on Nf2 and 4.8% on YM.

The gap is *not* launch overhead. CUDA API totals:

| API | Nf2 time (s) | Nf2 calls | YM time (s) | YM calls |
|---|---|---|---|---|
| `cudaMemcpyAsync` | 14.536 | 18612 | 2.445 | 3371 |
| `cudaDeviceSynchronize` | 0.675 | 427065 | 0.156 | 98227 |
| `cudaLaunchKernel` | 0.202 | 37099 | 0.040 | 6671 |
| `cudaMemcpy` | 0.120 | 33816 | 0.032 | 8652 |

`cudaMemcpyAsync` time (14.536 s / 2.445 s) tracks GPU kernel time (14.512 s / 2.434 s) almost
exactly — it is the blocking wait on the result copy, i.e. GPU time re-counted, not extra cost.
Everything else in the CUDA API is ~1.0 s (Nf2) / ~0.23 s (YM). So the missing ~7 s (Nf2) and
~2.4 s (YM) is **genuine host-side compute outside CUDA entirely** — the timestepper, spline
reconstruction, deal.II vector work. Note `cudaDeviceSynchronize` is called 427k times on Nf2
(598 per RHS evaluation) though it costs only 1.6 µs each.

## 4. Experiment 3 — LaunchBounds sweep {96, 128, 192}

Built with `-DCMAKE_CXX_FLAGS=-DDIFFRG_LAUNCH_BOUNDS=N -DCMAKE_CUDA_FLAGS=...` into
`build_lb96` / `build_lb192`; `build` is the stock 128. Verified before trusting anything:
the define reaches the ZA4 compile line and does **not** clobber DiFfRG's own flags (`-O3
-march=native -arch=sm_80` all still present), and the demangled kernel names carry
`LaunchBounds<96u>` / `<128u>` / `<192u>` — 44 kernels each.

### 4.1 Static: registers, spill, occupancy (sm_80)

Spill = STACK above the per-binary frame floor. That floor is code-specific — the small
`CT_get`-class kernels sit at 4000 B in Nf2 and 2416 B in YangMills — so it must be measured per
binary, not assumed. Occupancy is from registers at each variant's own block size.

QCD_Nf2:

| flow | LB96 REG/spill/occ | **LB128 (stock)** | LB192 REG/spill/occ |
|---|---|---|---|
| ZA4 | 168 / 2632 / 18.8% | **255 / 1864 / 12.5%** | 168 / 2632 / 18.8% |
| ZAAqbq2 | 168 / 1528 / 18.8% | **255 / 1040 / 12.5%** | 168 / 1528 / 18.8% |
| ZAAqbq1 | 168 / 1024 / 18.8% | **255 / 672 / 12.5%** | 168 / 1024 / 18.8% |
| ZAqbq4 | 168 / 2240 / 18.8% | **255 / 1680 / 12.5%** | 168 / 2240 / 18.8% |
| ZAqbq7 | 168 / 1968 / 18.8% | **255 / 1416 / 12.5%** | 168 / 1968 / 18.8% |
| ZAqbq1 | 168 / 1736 / 18.8% | **255 / 1168 / 12.5%** | 168 / 1736 / 18.8% |
| ZA3 | 255 / 56 / 9.4% | **255 / 56 / 12.5%** | 255 / 56 / 9.4% |
| ZAcbc | 254 / 48 / 9.4% | **254 / 48 / 12.5%** | 254 / 48 / 9.4% |
| ZA | 208 / 40 / 14.1% | **208 / 40 / 12.5%** | 208 / 40 / 9.4% |
| Zq / Zc | 168 / 16, 164 / 40 / 18.8% | **same / 18.8%** | same / 18.8% |

Two things stand out:

1. **LB=128 is the worst of the three for registers.** Both 96 and 192 push the six big kernels
   from 255 down to 168 regs, lifting occupancy 12.5% → 18.8%; 128 does not. 96 and 192 produce
   *identical* register and spill numbers — the shipped default sits in a bad spot between them.
2. The occupancy gain is bought with **~40–50% more spill** (ZA4 1864 → 2632 B). That is the whole
   trade, and only runtime settles it.

The small flows move the other way: ZA3/ZAcbc stay pinned at 255 regs and *lose* occupancy at both
96 and 192 (12.5% → 9.4%), because the odd block size wastes more of the register file. Together
they are 3.3% of GPU time, so this is second-order.

YangMills — the LB knob has almost nothing to work with:

| flow | LB96 | **LB128** | LB192 |
|---|---|---|---|
| ZA4 | 168 / 80 / 18.8% | **168 / 80 / 18.8%** | 168 / 80 / 18.8% |
| ZAcbc | 168 / 24 / 18.8% | **168 / 24 / 18.8%** | 168 / 24 / 18.8% |
| ZA3 | 128 / 24 / 23.4% | **128 / 24 / 25.0%** | 136 / 24 / 18.8% |
| ZA | 80 / 16 / 37.5% | **80 / 16 / 37.5%** | 80 / 16 / 37.5% |
| Zc | 72 / 16 / 42.2% | **72 / 16 / 43.8%** | 76 / 16 / 37.5% |

YangMills has **essentially no spill anywhere** (80 B on ZA4, i.e. one or two values) and its
dominant kernel is already at 168 regs / 18.8% at the stock 128. Nothing for LB to fix; the only
movement is ZA3/Zc getting slightly worse at 192.

### 4.2 Runtime — 4 alternated rounds, round 1 discarded, `fastInterpLookups` OFF throughout

QCD_Nf2:

| variant | median wall | vs LB128 | spread over kept rounds |
|---|---|---|---|
| **LB96** | **19.91 s** | **1.131x** | 0.07 s |
| LB128 (stock) | 22.52 s | 1.000x | 0.04 s |
| LB192 | 20.31 s | 1.109x | 0.04 s |

YangMills:

| variant | median wall | vs LB128 | spread |
|---|---|---|---|
| LB96 | 5.09 s | 0.982x | 0.02 s |
| LB128 (stock) | 5.00 s | 1.000x | 0.11 s |
| LB192 | 4.98 s | 1.004x | 0.04 s |

Physics gate: **identical to the last digit** across all six binaries
(Nf2 `1.1731490534726379` / `0.6888524203451547`, YM `6037.244194877355` / `0.9316318995037499`).
Step counts identical (350 / 240). Run-to-run spread is 0.02–0.11 s, so the Nf2 separation is
~60x the noise.

### 4.3 Result

**Nf2: LB=96 is worth 13.1% e2e** — by a wide margin the largest single effect in this session,
and 6.5x the 2% bar. LB=192 gets 10.9%. Since 96 and 192 produce *byte-identical* register and
spill numbers but differ by 2 points of wall, the win is not purely occupancy: at 128 threads the
5D `MDRangePolicy` tiling maps differently onto the 168-vs-255-reg allocation.

**YM: flat** (±1.8%, none clearing 2%) — as the static picture demanded, since YangMills has no
spill and is already at 18.8% occupancy at the default.

The important structural finding: **the shipped default of 128 is the worst of the three on
sm_80**, not a midpoint. That means the sweep straddled a bad point rather than bracketing an
optimum, so §6 extends it to LB ∈ {64, 160}.

## 5. Experiment 4 — `__noinline__` size-gate sweep {off, 300, 500, 800}

Variants built by rewriting `flows/*/kernels.hh` in place with `noinline.py`, building into
`build_ni<N>`, then `git checkout`-ing the headers back — rather than copying the source tree per
threshold. Verified per §8: the rewrite touches **only** decorator lines
(`KOKKOS_INLINE_FUNCTION` → `KOKKOS_FUNCTION __attribute__((noinline))`, 7 lines in nf2 ZA4 at
min=300, nothing else in the diff), the decorator counts match the static report at every
threshold, and each restore returned the tree to `noinline count = 0`.

⚠ The `git checkout` restore also reverts `flows/ZA4/numtrace.json`, undoing the `generated: 1`
pin from §1.1 — re-applied after each sweep, otherwise the next configure re-rolls the
nondeterministic `tr0` emit.

### 5.1 Static: the gate does **not** move registers on sm_80

| flow | off | ni300 | ni500 |
|---|---|---|---|
| ZA4 | 255 / 1864 B / 12.5% | 255 / **928 B** / 12.5% | 255 / 1896 B / 12.5% |
| ZAqbq4 | 255 / 1680 B / 12.5% | 255 / 1680 B / 12.5% | 255 / 1680 B / 12.5% |
| all others | unchanged | unchanged | unchanged |

This is the session's first clear refutation of a static prediction. The sm_90 cross-compile said
threshold 300 gives "0 spill / 168 regs / 18% occupancy". On **sm_80 none of that happens**:
registers stay pinned at 255 and occupancy at 12.5% for every threshold. The gate only *halves
ZA4's spill* (1864 → 928 B) at min=300, and 500 is a no-op (1896 ≈ 1864). Even ZAqbq4, which does
gate one function at 300 and 500, shows no spill change at all.

### 5.2 Runtime

QCD_Nf2 — flat:

| variant | median wall | vs off | spread |
|---|---|---|---|
| off | 22.39 s | 1.000x | 0.20 s |
| ni300 | 22.48 s | 0.996x | 0.08 s |
| ni500 | 22.48 s | 0.996x | 0.29 s |
| ni800 | 22.56 s | 0.992x | 0.08 s |

YangMills — clearly harmful:

| variant | median wall | vs off | spread |
|---|---|---|---|
| off | 5.00 s | 1.000x | 0.04 s |
| ni300 | 5.65 s | **0.885x** | 0.10 s |
| ni500 | 5.66 s | **0.883x** | 0.04 s |
| ni800 | 5.66 s | **0.883x** | 0.20 s |

Gate values pass everywhere. Nf2 `off`/`ni300` give `1.1731490534726379` and `ni500`/`ni800` give
`1.1731490534726354` — a last-ulp shift from different inlining, 2e-15 relative. YM is
bit-identical across all four.

### 5.3 Result — the headline experiment is a clean negative

| | 4070 | **A100** |
|---|---|---|
| Nf2 | 4.8% faster with the gate | **0.996x — nothing, ±0.4%** |
| YM | (not separately measured) | **0.885x — 13.0% slower with the gate** |

The gate should stay **off**, and the session doc's expectation ("300 best, 500 marginal, 800 ≈
off; if 300 wins clearly, add a 200 point") is not borne out — there is no win to extend, so no
200 point was built.

Why: out-of-lining pays only if it buys occupancy. On sm_90 it apparently does; on sm_80 ptxas
holds the big kernels at 255 regs regardless, so the gate buys nothing and charges call overhead.
YangMills shows the cost undiluted — its ZA4 has **no spill to recover** (80 B) and the gate
out-of-lines 65–82% of the flow's arithmetic, so it is pure loss. Nf2 breaks even only because
halving ZA4's spill roughly cancels the call overhead.

The follow-up worth testing is the one combination the sweep could not reach: the gate under
**LB=96**, where ZA4 spill is 40% higher (2632 B) — see §6.

## 6. Extended LaunchBounds sweep + the LB×noinline interaction

§4.3 showed the sweep straddled a bad point rather than bracketing an optimum, so three more
binaries were built: LB ∈ {64, 160}, and **LB96 + noinline 300** — the one combination the two
sweeps could not reach, and the only condition under which the gate could plausibly pay, since
LB=96 raises ZA4 spill 40% (1864 → 2632 B).

Static, QCD_Nf2 ZA4 (block size = each variant's own LB):

| LB | REG | spill | occupancy |
|---|---|---|---|
| 64 | 255 | 1864 B | 12.5% |
| **96** | **168** | **2632 B** | **18.8%** |
| 128 (stock) | 255 | 1864 B | 12.5% |
| 160 | 168 | 2632 B | 15.6% |
| 192 | 168 | 2632 B | 18.8% |
| 96 + ni300 | 168 | **1632 B** | 18.8% |

ptxas picks 255 regs at 64 and 128 but 168 at 96, 160 and 192 — the allocation is not monotonic
in the launch bound, and the two values that leave registers at 255 are exactly the two
powers of two straddling the winner.

Runtime, all six plus the stock binary alternated in one session, 4 rounds, round 1 discarded:

| variant | median wall | vs LB128 |
|---|---|---|
| LB64 | 20.84 s | 1.080x |
| **LB96** | **19.89 s** | **1.132x** |
| LB128 (stock) | 22.51 s | 1.000x |
| LB160 | 20.20 s | 1.114x |
| LB192 | 20.29 s | 1.109x |
| LB96 + ni300 | 19.91 s | 1.131x |

Gate values **bit-identical across all six**. Two conclusions:

1. **LB=96 is the optimum** and the curve is well-resolved around it (64: 1.080x, 96: 1.132x,
   160: 1.114x, 192: 1.109x). The stock 128 is a genuine local worst.
2. **The noinline gate is dead even where it should have helped.** Under LB=96 it cuts ZA4 spill
   by 38% (2632 → 1632 B) and changes wall by 0.1% (19.91 vs 19.89 s). Spill is simply not the
   bottleneck on this part at any launch bound — which retires the whole size-gate line of
   inquiry for sm_80, not just the {300,500,800} points.

## 7. Winner combination

**QCD_Nf2** — winner = `LB=96`, `fastInterpLookups` OFF, noinline off. Baseline = shipped state
(LB=128, flag defaults false). Both measured alternated in the §6 session:

| | wall | |
|---|---|---|
| shipped (LB128, flag OFF) | 22.51 s | 1.000x |
| **winner (LB96, flag OFF)** | **19.89 s** | **1.132x** |

**YangMills** — winner = `LB=128` (unchanged; 192 gains 0.4%, far under the 2% bar),
`fastInterpLookups` **OFF**, noinline off. Baseline = shipped state, which for YM is
`parameter.toml` with `fastInterpLookups = true`:

| | wall | |
|---|---|---|
| shipped (LB128, flag ON) | 5.04 s | 1.000x |
| **winner (LB128, flag OFF)** | **4.89 s** | **1.031x** |
| (LB192, flag OFF) | 4.87 s | 1.035x |

### 7.1 Winner profile

| | GPU kernel time | wall | GPU share |
|---|---|---|---|
| Nf2 baseline | 14.512 s | 22.51 s | 64.5% |
| **Nf2 winner** | **11.963 s** | **19.89 s** | **60.1%** |
| **YM winner** | **2.441 s** | **4.89 s** | **49.9%** |

The Nf2 win is entirely GPU-side: GPU time falls 2.55 s, wall falls 2.62 s, and the non-GPU
remainder is unchanged at 7.9 s. Per-flow, LB=96 helps exactly the kernels it was supposed to and
mildly hurts the mid-size ones:

| flow | baseline (s) | winner (s) | Δ |
|---|---|---|---|
| ZA4 | 4.748 | 3.665 | **−22.8%** |
| ZAAqbq2 | 4.040 | 3.064 | **−24.2%** |
| ZAAqbq1 | 3.182 | 2.492 | **−21.7%** |
| ZAqbq4 | 0.422 | 0.513 | +21.6% |
| ZAqbq7 | 0.389 | 0.419 | +7.7% |
| ZAqbq1 | 0.382 | 0.405 | +6.0% |
| ZA3 / ZAcbc | 0.257 / 0.223 | 0.256 / 0.226 | ~0 |

## 8. IBP go/no-go

⚠ **Note on the formula.** §9.2 of the session doc gives
`gain ≈ Σ_flow share(flow) × NTshare(flow) × 0.3` and quotes "~15% YM / ~35% Nf2 ceilings" from
the 4070. Those quoted numbers are reproduced by `Σ share × NTshare` **without** the 0.3
(4070 Nf2: 0.391×0.77 + 0.269×0.19 + 0.171×0.19 = 38.5%; YM: 0.869×0.20 = 17.4%). With the 0.3
the same inputs give 11.5% / 5.2%. So the quoted "ceiling" is the *NumTracer-arithmetic share of
GPU time*, and the 0.3 is the removable fraction applied on top. Both are reported below.
NTshares are the 4070 SASS attributions (nf2 ZA4 0.77, ZAAqbq1/2 0.16–0.22 → 0.19 used,
YM ZA4 0.20); they were not re-derived on A100.

| | NT-arith share of GPU time | GPU/wall | **IBP e2e ceiling** (0.3 removal) |
|---|---|---|---|
| Nf2 baseline | 34.6% | 64.5% | 6.7% |
| **Nf2 winner (LB96)** | **32.4%** (range 31–34%) | 60.1% | **5.8%** |
| **YM winner** | **16.6%** | 49.9% | **2.5%** |

**Verdict against the doc's stated criterion ("Nf2 ceiling ≥20% → GO"): GO — Nf2 is at 32.4%,
comfortably above 20%**, and barely moved from the 4070's ~35% despite fp64 going 1:64 → 1:2.
ZA4's share of GPU time fell (39.1% → 30.6%) but the ZAAqbq pair rose to compensate.

**But the honest e2e number is 5.8%, not 32%**, and it is falling for two independent reasons the
4070 could not show: the winner binary already removed 17.6% of GPU time, and only 60% of the Nf2
wall (50% of YM's) is GPU at all. An integration-by-parts reduction that removes 30% of all
NumTracer trace arithmetic buys **~5.8% e2e on Nf2 and ~2.5% on YangMills** on this hardware.
Recorded as instructed; **IBP not started**.

For comparison, the session's measured knob — a one-character change to `DIFFRG_LAUNCH_BOUNDS` —
is worth 13.2% on Nf2, i.e. **more than twice the entire optimistic IBP ceiling**.

## 9. Summary

| experiment | variant | Nf2 | YM | nf2 ZA4 regs/spill/occ |
|---|---|---|---|---|
| 1 | `fastInterpLookups` ON vs OFF | **0.953x** (ON costs 4.9%) | **0.986x** (ON costs 1.4%) | — |
| 3 | LB=96 vs 128 | **1.132x** | 0.982x | 168 / 2632 B / 18.8% |
| 3 | LB=192 vs 128 | 1.109x | 1.004x | 168 / 2632 B / 18.8% |
| 6 | LB=64 vs 128 | 1.080x | not built | 255 / 1864 B / 12.5% |
| 6 | LB=160 vs 128 | 1.114x | not built | 168 / 2632 B / 15.6% |
| 4 | noinline 300 | 0.996x | **0.885x** | 255 / 928 B / 12.5% |
| 4 | noinline 500 | 0.996x | 0.883x | 255 / 1896 B / 12.5% |
| 4 | noinline 800 | 0.992x | 0.883x | — |
| 6 | LB96 + noinline 300 | 1.131x | not built | 168 / 1632 B / 18.8% |
| 7 | **winner vs shipped** | **1.132x** | **1.031x** | |
| — | stock LB128 baseline | 22.51 s | 5.04 s | 255 / 1864 B / 12.5% |

### What the A100 changed, in one paragraph

Three of the four knobs the session was called to resolve came back **negative or reversed**, and
the one nobody flagged as interesting turned out to be the whole story. `fastInterpLookups`, worth
1.13x on the 4070, *costs* 1.4–4.9% here — the tabulated inverse transform is slower than the
transcendental it replaces once fp64 is 1:2, and it costs GPU time in every flow (+8.0%), worst in
the two interpolator-bound ones (+10%). The `__noinline__` size gate, billed as "the headline
experiment" and "the **only** lever that removes giant-kernel spill on datacenter parts", does not
move registers or occupancy at all on sm_80: it is neutral on Nf2, 13.0% slower on YangMills, and still
worth nothing under LB=96 where it does cut spill 38% (YM at 0.885x = 13.0% slower with the gate
on). What *does* work is the LaunchBounds value,
where the shipped default of 128 turns out to be a local worst on sm_80 — LB=96 is worth **13.2%
e2e on Nf2** with bit-identical physics. And the number that reframes everything downstream: with
kernels ~4.5x faster than on the 4070 and host work unchanged, **only 60% of the Nf2 wall and 50%
of the YangMills wall is still GPU**, so the next round's leverage is on the host side, not in the
kernels.

## 10. Decisions for the laptop (not made here — measurement only)

1. `DIFFRG_LAUNCH_BOUNDS` should become **arch-keyed**: 96 on sm_80, keep 128 on sm_89 (where it
   was tuned). This is the session's one large win and it is free.
2. `fastInterpLookups` default stays **false**; consider making YM's committed `parameter.toml`
   stop setting it true, or arch-key it — it is a 1.13x win on Ada and a 1.4% loss on Ampere.
3. The `__noinline__` size gate: **leave dead on sm_80.** The follow-on questions the doc parked
   as "premise void" (arch-keyed `gen.hpp` default, power table, caller-side trace-call hoisting)
   stay parked — the premise is void on this arch too, for a different reason (ptxas holds 255
   regs regardless, so out-of-lining buys no occupancy).
4. IBP is a **GO on the doc's stated criterion** but buys ~5.8% e2e on Nf2 / ~2.5% on YM. Worth
   weighing against host-side work, which is now 40% of the Nf2 wall and 50% of YM's and has never
   been profiled.
5. Not measured here, worth a future session: why `cudaDeviceSynchronize` is called 427k times in
   an Nf2 solve (598 per RHS evaluation), and what the ~7 s of non-CUDA host time actually is.

## 11. Artifacts

In `~/Code/numtracer/a100_logs/`:
`exp1_{nf2,ym}.tsv`, `exp3_{nf2,ym}.tsv`, `exp4_{nf2,ym}.tsv`, `exp5_nf2.tsv`, `win_ym.tsv`
(per-round wall/steps/gate values); `nf2_resusage_{base,lb64,lb96,lb160,lb192,ni300,ni500,lb96ni300}.txt`
and `ym_resusage_{base,lb96,lb192}.txt`; build/configure logs; `ZA4_regen.diff` and
`backup/nf2_ZA4_kernels.hh.orig` for §1.1.

nsys reports (in the run dirs, alongside their `*_cuda_gpu_kern_sum.csv` /
`*_cuda_api_sum.csv`):

| file | location |
|---|---|
| `nf2_a100_off.nsys-rep` | `QCD_Nf2/no_mesons/run_gate/` |
| `nf2_a100_winner_lb96.nsys-rep` | `QCD_Nf2/no_mesons/run_gate/` |
| `nf2_a100_on.nsys-rep` | `QCD_Nf2/no_mesons/run_fastinterp/` |
| `ym_a100_off.nsys-rep`, `ym_a100_winner.nsys-rep` | `YangMills/SP/run_gate/` |

Tools added to `gpubench-tools/`: `ab.sh` (alternating timing harness), `resusage.py` /
`rescmp.py` (register/spill/occupancy tables), `kernshare.py` (per-flow nsys shares),
`build_variants.sh` / `build_round2.sh` (the variant build drivers).

Variant build dirs left in place under each app (`build_lb*`, `build_ni*`, `build_lb96_ni300`) —
delete when done. Source trees are clean: `flows/**` restored to committed state
(nf2 ZA4 `kernels.hh` md5 `28d180cb57bf1c8bed8001f6db815605`), the only modified tracked file
being `flows/ZA4/numtrace.json` (`generated: 0` → `1`, the §1.1 pin).
