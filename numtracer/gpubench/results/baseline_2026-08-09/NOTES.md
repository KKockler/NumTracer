YM baseline 2026-08-09
qcd-codes @ dd46312 (pristine worktree), tune_STI flipped to false (only param change; matches production intent)
DiFfRG install: user rebuild 2026-08-09 19:55 from DiFfRG_KT 66a126e63
m2A = 6037.240528758876
Zc(0) = 0.9316141347304615
wall (single unwarmed run): 22.6 s

Nf2 baseline 2026-08-09
qcd-codes @ dd46312 (pristine worktree), tune_m2A flipped to false (only param change; m2A fixed at 8.767383098602 from config)
m2A = 1.17314905347264
Zc(0) = 0.688852420345153
wall (single unwarmed run): 106.5 s  <- first-ever Nf2 e2e reference

== WP1a result (Nf2 regen, 2026-08-09) ==
Regenerated all 11 flows with current emitter (interpshare+fma+cInline live; ZA4 kernels.hh 20573->10684 lines, 5688 fma).
Gate run (same fixed-m2A config as baseline):
  m2A = 1.1731490534726206 (rel diff 2e-14), Zc(0) = 0.6888524203451523 (1e-15)
  full output.h5: 1841 numeric datasets, worst max-rel-diff 3.07e-11  -> PASS (<1e-8)
  wall: 106.5 s -> 74.9 s (~1.42x, single unwarmed runs; protocol confirmation pending)
Note: production no_mesons/build was stale (cache pointed at renamed QCD_NumTracer path);
moved to build_stale_QCD_NumTracer_cache and reconfigured fresh.

== WP1b result (YM regen, Logarithmic, 2026-08-09) ==
Regenerated all 5 flows in the dd46312 worktree (current emitter; ZA4 kernel.hh: 10 shared _ix hoists).
Gate: m2A/Zc(0) identical to baseline at all printed digits; full output.h5: 1692 datasets,
420 bit-identical, worst max-rel-diff 7.2e-12 -> PASS.
wall: 22.6 s -> 17.4 s (single unwarmed runs).

== WP1c result (YM FocusedLog, production tree) ==
Regen+build+solve OK. m2A = 6037.244194877355 (6.1e-7 vs Log baseline), Zc(0) = 0.9316318995037505 (1.9e-5)
-> PASS (grid-truncation-level shifts, well under 1e-4 gate). Wall 20.0s (vs 17.4 Log; FocusedLog backward pricier -> WP2d target).

== WP1 e2e protocol (alternating, ratios only; absolutes drift with thermal state) ==
YM (Logarithmic, worktree): rounds 2-6 ratios 1.100/1.053/1.075/1.080/1.107 -> median ~1.08x
Nf2: rounds 1-2 ratios 1.063/1.082 -> ~1.07x  (the earlier single-run 1.42x was thermal-state confounded;
absolutes ranged 74.9-124s for the SAME binary+config across the day)
Conclusion: regen delivers ~1.07-1.08x e2e on both codes on the 4070 laptop; kernel-level gains are
diluted by CPU share of the solve. All physics gates PASSED (bit-level e-12/e-11 on fixed-grid regens).

== Post-regen SASS (probe layer) ==
results/l2_post_wp1_regen.csv: ym emission unchanged (50/265/1811/263/46 fp64) — nvcc had already
contracted; nf2 ZA4 probe 8453 fp64. Full-kernel attribution -> WP3.

== WP2a LANDED (2026-08-09) ==
Position cache in QuadratureIntegrator::map(): forward() once per grid point (device fill kernel,
same libdevice path) instead of per thread; key = to_string + size + hex first/last pos.
Gate: observables exact; 1688/1692 datasets bit-identical, worst 2.2e-15 (lost cross-inline fma
contraction of forward into kernel — benign). SASS per map kernel: 4440->4420 / 2571->2546 /
1864->1839 fp64-core. Wall: neutral within 4070 noise (expected; payoff on A100 + FocusedLog).
NVCC LESSON: if-constexpr INSIDE extended lambdas in this templated member = miscompiled
(cudaErrorInvalidResourceHandle or silent zero RHS); lambda DEFINED inside if-constexpr = fine.
Pattern used: two complete branch-free lambdas under host-level if constexpr.

== WP2b LANDED (2026-08-09) ==
k-only lookup hoisting: Codegen.m expression-level pass (HoistLoopConstLookups, on via DiFfRG_compat
when interp types match; NT_GEN_NO_KHOIST=1 A/B) + ntHoisted static host evaluator + compat patches
CT_map_*/CT_get wrappers AND the tuple-forwarder (by-value unpack broke .CPU() on copies).
YM: 6 lookups hoisted per flow; SASS map kernels -451..-455 fp64 each (ZA4 -10.2
== WP2b LANDED (2026-08-09) ==
k-only lookup hoisting: Codegen.m expression-level pass (HoistLoopConstLookups, on via DiFfRG_compat
when interp types match; NT_GEN_NO_KHOIST=1 A/B) + ntHoisted static host evaluator + compat patches
CT_map_*/CT_get wrappers AND the tuple-forwarder (by-value unpack broke .CPU() on copies).
YM: 6 lookups hoisted per flow; SASS map kernels -451..-455 fp64 each (ZA4 -10.2%, ZA3/ZAcbc -17.7%,
ZA/Zc -24.7%). Physics: worst rel 5.3e-11 over 1692 datasets (host-vs-device last ulp). Wall on 4070:
+1.5-2.5% e2e (r2 1.015, r3 1.023) - A100 upside is the real target. ctest 41/41 green.
Debug lessons -> memory: NumTracer` context Protected (Set::wrsym silent); scaffold tuple-forwarder
unpacks by value.

== WP2b production landing ==
Nf2: 11 flows hoisted; gate solve m2A = 1.1731490534726272 (6e-15 vs baseline), Zc(0) 4e-16. PASS.
YM (FocusedLog): 5 flows x 6 lookups; m2A = 6037.244194877355 IDENTICAL to WP1c ref, Zc(0) last ulp. PASS.

== WP2d LANDED (2026-08-10) ==
Tabulated inverse index() behind runtime ConfigTree flag /discretization/fastInterpLookups
(default false; sampled once at interpolator construction via ConfigurationHelper).
Implementation: per-cell Chebyshev fit of backward() (adaptive degree 6..16, target 1e-12,
throw >1e-10), Clenshaw eval + 6-step branchless binary search (~50 device ops vs ~155
log1p / ~350-400 log+asinh). Standalone math check: deg 8 @ 2-4e-13 (Logarithmic), deg 14
@ 2e-13 (FocusedLog). Unit test [fastindex] added to DiFfRG spline_interpolator_1d.cc:
8 assertions pass; full spline suite 36863 assertions green.
YM production (FocusedLog): flag-off IDENTICAL to WP2b gate; flag-on m2A rel 6.3e-15.
WALL: 1.134x / 1.129x (stable, off 17.4-17.6s vs on 15.3-15.4s) -> ~1.13x e2e, the largest
single lever this round. Nf2: gates pass (8.8e-14); wall pending.
TRAP (recorded in memory): cmake --install preserves source mtimes; consumers silently
skip rebuilding — touch installed headers + verify "Building CXX" count.

== WP3 attribution (attrprobe.py, full CT_map TU, cumulative stubbing; static fp64, sm_89, map kernel fn) ==
flow        full   NT-in-kernel  interp  regulators  remainder
ym ZA4      4830   1833 (38%)    2181*   516 (11%)   300
nf2 ZA4    14212  10969 (77%)    2141*   567 (4%)    535
nf2 ZAqbq4  7463   3561 (48%)    2022*   573 (8%)   1307
(* interp bucket includes BOTH index paths — the runtime fastInterpLookups flag keeps analytic+table
   compiled; runtime cost is one path.)
The former "unattributed 50-68%" = interpolator at()+transforms + regulators(4-11%) + small remainder.
=> WP2c (l1-only per-node precompute) REJECTED: its slice is ~1/3 of the regulator bucket, ~1% wall.
=> nf2 ZA4 map kernel is 77% emission arithmetic: WP4 (power table/Horner) and any future IBP aim there.

== WP3: first Nf2 GPU profile (nsys over the fixed-m2A solve, 714 RHS) ==
ZA4 39.1% | ZAAqbq2 26.9% | ZAAqbq1 17.1% | ZAqbq4/1/7 ~3% each | ZA3 1.6% | rest <1.5%
=> ZA4 + ZAAqbq1/2 = 83.1% of GPU time. The ZAqbq family (9.3%) is NOT the target the old
   65%-figure suggested; ZAAqbq1/2 (44%) were invisible in the old accounting.
eta loop: 6.34 iterations/RHS (4524/714) — ZA/Zc/Zq still only ~4% total.
Attribution of the two: ZAAqbq2 = 57% interp / 22% NT / 12% reg; ZAAqbq1 = 58% interp / 16% NT.
=> interp-dominated: fastInterpLookups (WP2d) is their lever (and already landed);
   nf2 ZA4 (77% NT) is the WP4 emission target.
Full data: results/attr_2026-08-10.csv, run_gate/nf2_profile.nsys-rep.

== WP4 CLOSED (2026-08-10) ==
- Horner orders: NT_GEN_HORNER_ORDERS=8 on nf2 ZA4 -> 10684->10670 emitted lines (-0.13%).
  REJECTED (gate >=2%). The gen.hpp "~1%" comment holds on the giants too.
- Power table: REJECTED WITHOUT IMPLEMENTATION. Premise (cross-noinline-boundary redundancy)
  is void: grep noinline = 0 in fresh production kernels (the ntKokkosDecor decorator trap
  keeps the gate inert in production), so all traces inline and nvcc GVN already shares the
  power products. COUPLING: if the noinline gate is revived for the sm_90 spill fix (WP5),
  re-evaluate the power table.
- IBP: remains open, now with real weights: nf2 ZA4 = 39% GPU x 77% NT-in-kernel; ym ZA4 =
  87% x 38% (static shares, dual-index-path inflated). Decide after A100 runtime.
- SIDE FINDING: kernels.hh bytes depend on the generator WORKER COUNT (term-grouping order):
  numtrace under -j8 x 11 flows vs solo gives different-but-equivalent tables (verified: gate
  solves agree to 2e-16; deterministic at fixed worker count). Never byte-diff kernels.hh
  across builds with different parallelism.
