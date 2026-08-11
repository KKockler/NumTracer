# A100 session 2 — handoff

Successor to `NUMTRACER_A100_SESSION.md` → `A100_NOTES.md`. Same discipline: alternate variants,
discard round 1, gate values every run, ratios not absolutes.

**Read first:** `docs/NUMTRACER_TILE_WASTE.md` and `docs/NUMTRACER_PER_THREAD_FRAME.md`. Between
them they retire the headline conclusion of session 1.

---

## 0. What changed since session 1, and what it invalidates

Two defects were found and fixed in `DiFfRG_KT` on the laptop. Both are in the launch path, neither
is in NumTracer.

**(a) The MDRange tile masked a third of the threads.** `Kokkos::LaunchBounds<maxT>` clamps the tile
as well as hinting registers, and Kokkos can only halve its recommended tile — so every tile
dimension it produces is a power of two. Every QCD_Nf2 model uses angular order **6**, which no
power-of-two tile divides. At the shipped LB=128 the rank-5 tile `{16,2,4,1,1}` launches 8 slots for
6 points, and since dim0×dim1 = 32 is one warp, those are whole dead warps.

**This is what session 1's LaunchBounds sweep actually measured.** Cross-checked against §7.1
per-flow: the three rank-5 flows (the only ones with waste) gained 21.7–24.2% against 25.0%
predicted from the masked fraction alone; every rank-4 flow — same occupancy change, nothing to
recover — came out **6–21.6% slower**. Clinchers: LB=64 beat LB=128 by 8.0% at *byte-identical*
registers, spill and occupancy; and YangMills, the one code with angular order 8, was flat across
the whole sweep.

> **Therefore: do not act on `A100_NOTES.md` §10.1 ("arch-key `DIFFRG_LAUNCH_BOUNDS` to 96 on
> sm_80"). That recommendation is withdrawn.** With the tile fixed, LB=128 may be correct
> everywhere. Re-measuring it is item 4 below.

Masked fraction at LB=128: YangMills **1.000x**; nf2 no_mesons 1.255x; with_mesons 1.306x;
with_mesons_3D 1.299x; finite_T no_mesons 1.314x. Recompute for any app with
`numtracer/gpubench/tools/tilewaste.py <app_dir>` — it needs no GPU.

**(b) Every thread carried a by-value copy of all 13 interpolators.** `device::tuple_cat(x, pos,
m_args)` built a fresh per-thread tuple in local memory. That is the "per-binary frame floor" §4.1
mentions and dismisses: **3576 B of the 4000 B frame**, twice the size of the spill two rounds went
into. Fixed by passing nested packs (references) instead.

Measured on the **4070**: **2.323x**. The frame was being paid for, not merely allocated. Untested
on sm_80 — that is item 2.

**The two fixes have opposite architecture dependence, and that is the point:** the tile fix is
1.001x on sm_89 and (inferred) 23% on sm_80; the frame fix is 2.3x on sm_89 and unknown on sm_80.
Masked warps only cost you when the SM is short of eligible warps; a local-memory frame costs you
whenever you touch it. Do not assume either transfers.

---

## 1. Build state

| | |
|---|---|
| DiFfRG | `DiFfRG_KT` @ `feature/fast_interpolators`, HEAD `d8e10c83a` + the two fixes |
| relation to session 1 | session 1 built `7ab9cc5fa`; `d8e10c83a` ("async Maps") is its child |
| does that matter? | `async Maps` did **not** touch the device launch path — only the host-side async result copy. It *did* remove `fastInterpLookups` entirely (deleted `fast_interp.hh`). Since session 1's verdict was "flag OFF", HEAD **is** the measured flag-OFF baseline. |
| switches | `DIFFRG_DIVISIBLE_TILE` (default 1, `=0` control). The frame fix is **unconditional** — no flag. |

Verify before trusting anything: `numtracer/gpubench/tests/tile/run.sh` must pass, and
`nvcc --version` on the cluster is 12.8 vs 12.9 locally — reproduce **nf2 ZA4 → 255 regs / 1864 B at
LB=128, sm_80** once before believing any static number carried over.

---

## 2. Agenda

Ordered so the long build runs in the background and the decisive measurement comes early.

### 1. Start the `lambda1L3D` `-Xptxas -O3` build (background, hours)

`with_mesons_3D/CMakeLists.txt:9-28` pins `lambda1L3D`'s device TUs to `-Xptxas -O1` **only because a
30 GB laptop OOM-killed the `-O3` build at ~22 GB**. The node has 503 GB. Remove the pin, build,
record ptxas peak RSS / wall / registers / spill, and compare runtime against the `-O1` build. No
patch needed. This may be worth more on that code than every occupancy knob combined.

### 2. `ncu --set full` — the measurement no static work substitutes  ⟵ **highest value**

Nsight Compute has **never** been run on these kernels; every prior round used `nsys` (timeline) and
static `cuobjdump`. Targets: nf2 ZA4, ZAAqbq1, ZAAqbq2, YM ZA4.

Take it on the **fixed** build, and on a build with the frame fix locally reverted (revert the three
`device::apply` sites in `quadrature_integrator.hh` back to `tuple_cat`). What to read:

- `sm__pipe_fp64_cycles_active` — settles whether these kernels are fp64-pipe-bound (occupancy is a
  red herring) or latency-bound (it isn't). The laptop roofline estimate suggests nf2's big three sit
  near ~50% of fp64 peak and YM ZA4 near ~90%, which would retrodict both session-1 outcomes; confirm
  or kill it.
- local-memory throughput and `long_scoreboard` stalls, before vs after the frame fix — how much of
  the 2.3x survives at fp64 1:2.

**If fp64 utilisation comes back >70%:** occupancy is not the lever, items 5–6 are pointless, and the
next move is arithmetic reduction (IBP) or the host side. Learning that in 15 minutes is a good
outcome, not a wasted session.

### 3. Frame fix, end-to-end

Fixed vs locally-reverted, both apps, alternated 4 rounds. This is the one whose sm_80 magnitude is
completely unknown — 2.3x on the 4070 is an upper bound, since that probe has far less arithmetic to
hide the traffic behind than a real flow.

### 4. Re-measure LaunchBounds **with the tile fixed**

`DIFFRG_LAUNCH_BOUNDS ∈ {64, 96, 128, 192}` × `DIFFRG_DIVISIBLE_TILE=1`. Session 1's sweep is void:
it was measuring divisibility. Expect the LB curve to be much flatter now, and expect 128 to stop
being a local worst. Also run `DIFFRG_DIVISIBLE_TILE=0/1` at LB=128 directly — the clean A/B for the
tile fix itself, predicted ~1.20x on nf2 GPU time and ~1.0x on YangMills.

### 5. `minBlocksPerMultiprocessor` — now a *low-prior* experiment

`Kokkos::LaunchBounds<maxT, minB>` — the second parameter is still unset, and it is the only way to
move registers without touching the tile (`-maxrregcount` is inert; `__launch_bounds__` outranks it).
But the rank-4 control in §0 says raising occupancy 12.5→18.8% was a *net loss* when it wasn't
bundled with a tile fix, so the prior is much weaker than session 1 suggested.

Only run this if item 2 says the kernels are latency-bound. If so: hold `maxT=128`, sweep
`minB ∈ {3,4,5}` (reg caps 168/128/102, i.e. 18.8%/25%/31%). Apply `minB` to the `dim>1` MDRange
specialization **only** — the `dim==1` `RangePolicy` feeds `cuda_deduce_block_size`, where `minB` can
drive the block size to 0. Key any default on `KOKKOS_ARCH_*`, never `__CUDA_ARCH__` (the tile clamp
is host-side; an arch-dependent policy type is an ODR hazard that compiles cleanly).

### 6. `fastInterpLookups` × launch bounds — never tested together

Session 1 measured the flag at LB=128 only and LB at flag-OFF only. The flag lost 10% *at 12.5%
occupancy*; if its cost was unhidden memory traffic it may flip at higher occupancy. Note it no
longer exists at HEAD (removed by `async Maps`), so this needs `7ab9cc5fa`. Lowest priority.

### 7. Winner combination, then `nsys` for the per-flow table

Mirror `A100_NOTES.md` §7.1 so the two sessions are directly comparable.

---

## 3. Gates and traps

Gate values — **check every run**, they have been bit-stable across every variant so far:

| code | m2A | Zc(0) | steps |
|---|---|---|---|
| QCD_Nf2 | 1.1731490534726379 | 0.6888524203451547 | 350 |
| YangMills | 6037.244194877355 | 0.9316318995037499 | 240 / 425 RHS |

Expected numeric behaviour of the new changes: the **tile** change is bit-identical (a `parallel_for`
with one write per thread — re-tiling reorders no arithmetic; verified over 64 points). The **frame**
change is *not* — it shifts results by ~1 ulp (1.3e-16 relative), because re-binding arguments
changes what nvcc inlines and hence how it contracts FMAs. Same class as the `ni500`/`ni800` shift in
§5.2. Tolerance 1e-8 as before.

Traps, all of which cost time in session 1:

- `git checkout` of the flow headers **reverts the `flows/ZA4/numtrace.json` `"generated": 1` pin**;
  re-apply after every restore or the next configure re-rolls a nondeterministic `tr0` emit.
- `cmake --install` preserves source mtimes → consumers silently skip rebuilding. `touch` the
  installed headers and confirm "Building CXX" actually fires.
- The GPU idles at 210 MHz; whichever variant runs first looks 20–40% slow. Alternate, discard
  round 1, report ratios. `gpubench-tools/ab.sh`.
- `ctest -L codegen` passes green even when generation aborted.
- CUDA 12.8 (cluster) vs 12.9 (laptop): reproduce one known register/spill number before trusting
  any static figure carried over.

---

## 4. Do not queue

`-maxrregcount`; the `__noinline__` size gate at any threshold; `fill()` CSE; `CrossTraceCSE`;
CUDA-graph batching of the eta loop; explicit tiles matched to the launch-bounds *product* (the
divisibility-aware tile is a different thing and is already in). All measured, all inert or harmful.

## 5. Still unexamined, and probably the largest remaining block

**Host-side work has never been profiled.** GPU is 64.5% of the Nf2 wall and 47.9% of YM's; the
remainder is ~7.9 s / ~2.4 s of genuine host compute — not launch overhead (`cudaLaunchKernel` totals
0.202 s over 37k calls). It includes 427,065 `cudaDeviceSynchronize` calls per Nf2 solve, 598 per RHS
evaluation against ~52 kernel launches per RHS, which smells like a synchronising accessor inside the
eta fixed-point loop.

This needs **no cluster time** — host code is architecture-independent, so `perf record -g` on the
laptop profiles the same timestepper at the same problem size. It is the one measurement in this
programme that transfers at full fidelity, and it competes with nothing on the agenda above.
