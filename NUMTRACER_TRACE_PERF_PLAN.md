# Large-trace generator efficiency — staged implementation plan

## North star: shrink RAM so ZAAqbq can run at a HIGH phase-B thread count

The point of the RAM work is not just to fit in memory — it is to **let the ZAAqbq flows run with many
phase-B workers (`nB`) instead of being forced down to `nB=1..2` to avoid OOM.** Phase B does the
CSE/Horner *lowering*, which is compute-heavy and parallel, so a low `nB` cripples runtime: measured on
the pure-gauge `za4_147`, phase-B+lower is **625 s at `nB=2` vs tens of seconds unthrottled**. Each
phase-B worker holds one full contraction, so today ZAAqbq must throttle `nB` low for RAM — which is
exactly what makes it slow. Every term-size / peak reduction below (1c −12.5%, Stage 2 −44%, Stage 3
blocked `operator*`) directly raises the `nB` that fits in a given RAM budget, and thus the speed. RAM
is the means; **a high-`nB` ZAAqbq run is the goal.**



Companion to `NUMTRACER_TRACE_PERF_FINDINGS.md`. Every step is landed against a **small** benchmark
(never the 12-min full ZAAqbq regen in the inner loop), quantified before/after, and correctness-gated
before moving on. The full flow is regenerated only once per batch, to confirm the headline.

## Measurement ladder (cheapest → most expensive)

| L | Harness | Cost | Measures | Correctness |
|---|---|---|---|---|
| **L0** | Op-level micro-bench against the real headers (the investigation's `bench_mul.cpp` / `pb_bench.cpp` / `bench.cpp` / `probe.cpp`, kept in the scratchpad) | seconds | the isolated op's RAM peak / time / `sizeof` — this is where the *magnitude* of a RAM lever shows (the small flow is too small to reach the >8 GB collapsing regime) | `eval()`/checksum bit-identity built into each bench |
| **L1** | `cmake --build build && ctest --test-dir build` (no Wolfram) | ~1–2 min | — | the full default suite must stay green (`test_trace_fold`, `test_numeric_contract`, `test_dpoly`, the `compare_*` flows against committed headers) |
| **L2** | **Standalone small generator**: recompile the committed `tests/gen/gen_zaaqbq1_small_num.cpp` + `_u*.cpp` against the patched headers, run under `NT_GEN_PROFILE=2` | ~1–2 min | end-to-end generator RAM + phase A/B time on a real dense-flow proxy (7,500 contractions / 820 traces / longest fold 140) | `diff` emitted header vs committed `ZAAqbq1Small_num_kernel.hh` (expect bit-identical for pure-refactor levers), then `compare_zaaqbq1_small` for the scale-relative grade |
| **L3** | `ctest -L codegen` (regenerates every flow from Wolfram + grades) and/or full ZAAqbq1 standalone gen under a `systemd-run … MemoryMax` cap | ~7–12 min | the **headline** RAM/time on the true monster | the codegen grade + `MemoryPeak`. Run **once per landed batch**, not per iteration |

### L2 recipe (the workhorse — no Wolfram)
```sh
cd numtracer
g++ -std=c++20 -O1 -w -Iinclude -Itests/gen \
    tests/gen/gen_zaaqbq1_small_num.cpp tests/gen/gen_zaaqbq1_small_num_u*.cpp -o /tmp/gen_small
# read the committed -n/-d off the .wls or the header's namespace; then:
NT_GEN_PROFILE=2 NT_GEN_MAXW=6 NT_GEN_MAXW_B=1 /tmp/gen_small -n <ns> -d <decor> > /tmp/small.hh
diff <(sed -n 's/.*//p' /tmp/small.hh) tests/gen/ZAAqbq1Small_num_kernel.hh   # bit-identity gate
```
The generator's own `[cse]`/profile lines give the distinct-trace count, phase-A/phase-B seconds, and
trace-table MB. For L0 RAM peaks use `/usr/bin/time -v` (Max RSS) or a `systemd-run --user --scope
-p MemoryMax=…` wrapper; per `NUMTRACER_LOWERING_RAM.md`, trust systemd `MemoryPeak`, not the per-wave
RSS lines (they miss mid-wave spikes).

### Non-negotiable gotchas (from `NUMTRACER_LOWERING_RAM.md`)
- **Rebuild `libNumTracer.a` before any L3 Wolfram regen** — `trace_fold.hpp` is header-only and not in
  the archive; `touch src/numtracer.cpp` first or `MakeNTKernel::stalelib` aborts (or silently links a
  stale vintage → wrong traces in a byte-identical-looking header).
- **A codegen pass can be false-green.** Check `genfail == 0` and `wrote generator` in each
  `/tmp/regen_gen_*.log`, not the ctest exit code.
- **Grade on numbers, not bytes.** The dedup can renumber `sN` while every value is identical; and the
  kernel crosses zero so pointwise *relative* error is meaningless near a root. `compare_zaaqbq1_small`
  already grades on scale-relative `max|a−b|/max|k|` (target < a few eps) + median pointwise deviation.

## Stage 1 — the "free" batch (near-zero risk, land together)

Order chosen so each change is independently revertible and independently measured.

### Step 1a — `MPolyFactory::scaled` (lever B1)
- **Change.** Add a factory that scales every coefficient of an `MPoly` by a constant `Cx` in place,
  preserving term order (no scratch, no `std::sort`). Repoint the three call sites: `scale_trace`
  (`trace_fold.hpp:61`), the colour-fold scale, and the `DPoly` `scaleCx` path (`dpoly.hpp:167`).
- **Quantify (L0).** `pb_bench.cpp` — expect 2.8–6.3× on the scale op (DPoly 6-dmono×60-term:
  ~745→~157 ms). **Correctness:** bench checksum bit-identical.
- **Gate (L1 + L2).** `ctest` green; L2 emitted header **byte-identical** to committed (this is a pure
  refactor — any diff is a bug). Record phase-B seconds from the L2 profile line; expect ~15–30% lower.
- **Rollback if:** header not bit-identical, or phase-B time not improved.

### Step 1b — `tcmalloc_minimal` preload (lever T)
- **Change.** Link the disposable generator binary with `-ltcmalloc_minimal` (via the codegen g++
  command in `Codegen.m`, gated on availability), or document `LD_PRELOAD` for the standalone path.
  No source logic changes.
- **Quantify (L2).** Same `/tmp/gen_small` run with/without preload — expect ~−7% wall time.
  **Correctness:** header byte-identical (allocator swap cannot change output).
- **Rollback if:** not available or no measurable gain (keep it opt-in, never a hard dep).

### Step 1c — 32-bit `size_type` on the `Mono` small_vectors (lever D2)
- **Change.** Give `MonoExp`/`MonoAtoms` a `size_type`-narrowed allocator/traits (24→16 B header);
  optionally tighten the atom inline count. Two `using` lines, no arithmetic touched.
- **Quantify (L0).** `probe.cpp` — `sizeof(pair<Mono,Cx>)` 128→112 (−12.5%); feed through `poly_bytes`
  to project the net-poly and single-contraction floors. **Correctness:** bench checksum bit-identical.
- **Gate (L1 + L2).** `ctest` green; L2 header byte-identical; note L2 peak RSS delta.
- **Rollback if:** any value changes, or a small_vector API assumption breaks.

**Batch exit (one L3).** After 1a–1c: `touch src/numtracer.cpp`, rebuild lib, `ctest -L codegen`
(all green, `genfail==0`), then one full ZAAqbq1 standalone gen under `MemoryMax` for the headline
RAM/time. Record in the findings doc.

## Stage 2 — packed-exponent `Mono` (levers A1 + D1, the structural win)

Do only after Stage 1 is green — it is the one change that touches every `e[k]` site.

- **Change.** Pack exponents into `u64[2]` (24 syms × 5 bits): (i) storage shrinks the term (128→72 B,
  or 56 with tightened atoms); (ii) `from_scratch` sorts on the packed `u64` key instead of walking the
  exponent small_vector in `Mono::operator<`. **Keep an exact element-wise `e`-fallback** for
  degree>31 or nsym>12 (both agents insist on this — it makes the change exact under overflow).
- **Quantify.** L0 `probe.cpp` for the new `sizeof`; L0 `bench.cpp` for the sort speedup (expect
  −7…9% per contraction, and the collapsing-regime peak in `bench_mul.cpp` down proportionally to the
  term shrink, >8 GB → ~4.5 GB projected). **Correctness:** every bench checksum bit-identical.
- **Gate.** L1 `ctest` green (esp. `test_numeric_contract`, `test_dpoly` — these exercise Mono
  ordering directly); L2 header **byte-identical** (packing must not reorder like terms); L3 codegen
  grade + full-flow RSS.
- **Rollback if:** any header diff, any bench value change, or the overflow fallback is ever wrong on a
  high-degree probe (add a targeted degree-32 unit case first).

## Stage 3 — blocked `operator*` (lever R1, caps the last transient)

- **Change.** In `MPoly::operator*` (`mpoly.hpp:262`), gate on scratch size: below `kMulMaxScratch`
  (~2²⁰) the current path is untouched (protects the mean-4-term hot path); above it, block the outer
  operand so live scratch is capped at `rows·|b|`, blocks merged via `operator+`.
- **Quantify (L0).** `bench_mul.cpp` collapsing regime — the harness already shows the 1944 MB peak;
  expect it capped to a few hundred MB (15–25× on the pathological multiply). **Correctness:** checksum
  bit-identical; the only change is a block-sum reassociation (< the shipped tree fold's ~7e-15).
- **Gate.** L1 `ctest` green; L2 header — expect byte-identical (small flow stays under the threshold,
  so the gate proves the fast path is untouched); a targeted L0 above-threshold case proves the blocked
  path matches the unblocked value to tolerance.
- **Rollback if:** hot-path (below-threshold) result not bit-identical, or the block path exceeds the
  reassociation tolerance.

## Stage 4 — eager quark-gluon vertex collection (lever W1, the root-cause fix)

Separate, larger track — retires `dropAqbq47` with the full ZAqbq1/4/7 basis. Two sub-steps, the first
of which already wins the generator-RAM/size metric on its own.

- **4a (emitted-size, exact).** Extend `collectibleDiracSumQ`/`diracNumeratorSumQ` (`DSL.m:221-234`) so
  a multi-structure vertex with an open Lorentz leg is collected (`ntVertexSlot` carrying the free
  gluon index) instead of distributed. **Quantify (L2 + L3):** distinct-trace count and table MB on the
  *small* flow first, then ZAAqbq1 (expect 488→~24 nets, 20.1→~1–2 MB). **Correctness:** this changes
  the emitted set, so grade with `compare_zaaqbq1_small` (scale-relative), **not** byte-identity — and
  cross-check the collected result against the current distributed `dropAqbq47`-free path on a flow
  small enough to generate both.
- **4b (contraction-RAM follow-on).** Carry the vertex slot through contraction as a matrix-of-`DPoly`
  (the `DSlotOpt` "later extension", `numeric_contract.hpp:943-951`) so the eager form also shrinks the
  contraction working set, not just the emitted set. Larger; scope after 4a is measured.
- **Rollback if:** the collected kernel disagrees with the distributed reference beyond grading
  tolerance on any small flow.

## Stage 1 — MEASURED RESULTS (landed 2026-07-22)

Steps 1a and 1c are landed in `include/numtracer/numeric/{mpoly,dpoly,trace_fold}.hpp`. 1b (tcmalloc)
is available but not wired in (see below).

### 1a — `MPolyFactory::scaled` (kill the dead sort in constant scaling) — LANDED
- **L0 correctness:** bit-identical to `constant(ns,c) * p` across every term count + the `c==0` edge,
  and DPoly `scaleCx` across every case (`scratchpad/l0_scaled.cpp`, memcmp on the coeff bits).
- **L0 speed:** MPoly scale **3.7–4.0×** (154→41 ns at 4 terms, 1159→291 ns at 60); DPoly `scaleCx`
  **4.1×** (7105→1717 ns on the 6-dmono×60-term case). Squarely in the predicted 2.8–6.3× band.
- **End-to-end neutrality (control):** regenerated ZA3 trace table WITH vs WITHOUT the edit is
  **byte-identical**; the ZAAqbq1Small kernel regenerated `unchanged`. So the emitted kernels do not
  move — the win is pure generator-time.

### 1c — 32-bit `size_type` on the `Mono` small_vectors (`NarrowAlloc`) — LANDED
- **RAM (probe, `scratchpad/probe.cpp`):** `MonoExp` 72→64, `MonoAtoms` 40→32, **`Mono` 112→96**,
  **stored term `pair<Mono,Cx>` 128→112 (−12.5%)** — so every `poly_bytes` floor drops 12.5% (the
  ~41 MB net poly → ~36 MB; the >8 GB single contraction → ~7 GB, which may itself clear the 10 GB cap).
- **Correctness:** inline/heap thresholds unchanged (`resize(24)` inline, `resize(30)` heap — only the
  bookkeeping width shrank, never a value); the whole engine + all tests recompile; **all 35
  non-codegen tests pass** (`ctest -LE codegen`), including the `compare_*` value-grading flows, so no
  computed coefficient moved; and 1a's L0 bit-identity still holds with the narrowed `Mono`.

### 1b — tcmalloc preload — AVAILABLE, not wired
`/usr/lib/libtcmalloc_minimal.so` is present. Not added to the generator compile because that lives in
`Codegen.m`, which currently carries unrelated in-flight edits; realise it non-invasively with
`LD_PRELOAD=/usr/lib/libtcmalloc_minimal.so` on the generation run (agent-measured ~−7%). Wire it into
`Codegen.m`'s g++ command (gated on availability) once that file settles.

### Note on the `codegen_regen` ctest
`codegen_regen` is **pre-existing red**, independent of Stage 1: the committed `tests/gen/*.hh`
fixtures are stale against the working tree's modified `Codegen.m`/`DiFfRG_compat.m` (regenerating any
flow, e.g. ZA3, changes them the SAME way with or without my edits — cleaner FP constants like
`13.499999999999943 → 13.5`), and the density guard reads a `[cse]` line the NDIAG-cut zaaqbq flow no
longer emits. The clean end-to-end byte-identity + headline phase-B/RAM (L3) should be re-measured once
the fixtures are refreshed against a settled `Codegen.m`.

## Stage 2 — MEASURED RESULTS (packed-exponent Mono, landed 2026-07-22)

`MonoExp` is now a **small-buffer-optimized packed key**: 2 inline 64-bit words hold up to 24 symbols at
5 bits each (big-endian per word, so the integer order of `w[0]` then `w[1]` reproduces the OLD
element-wise lexicographic exponent order EXACTLY), with a heap `std::vector<MonoExpT>` fallback the
moment a symbol index ≥ 24 or an exponent > 31 appears. So nsym and degree are **UNBOUNDED** (whatever
flow the user provides), exactly like the old heap-expanding small_vector — NOT a hard cap.

- **RAM:** `Mono` 112→**56 B**, stored term `pair<Mono,Cx>` 128→**72 B (−43.75%)** — smaller than a
  fixed 4-word cap would have been (80 B), because 2 words + one pointer (24 B) beats 4 words (32 B).
  The >8 GB single contraction → ~4.5 GB; 41 MB net polys → ~23 MB. This is the term shrink that raises
  the `nB` fitting a RAM budget — the north star.
- **Runtime:** on a contraction microbench (operator* + from_scratch sort, nsym=12, interleaved,
  `-march=native`) Option B is **~42 µs/mul vs the original small_vector's ~61 µs — ≈30% FASTER**, with
  identical results (`sink` matches bit-for-bit). The non-trivial copyability does not eat the win: the
  common-case copy/move is 3 words + a predicted null-check, while the hot `operator<` drops from an
  element walk to a single word compare.
- **Correctness:** round-trip get/set OK and packed ordering IDENTICAL to lexicographic across inline
  AND heap (720k pairs, nsym up to 48); all 35 non-codegen tests pass INCLUDING the numeric_contract
  stress test tr(p̸₁…p̸₁₀) at nsym=40 which now exercises the heap fallback; 1a's scaled bit-identity
  still holds. End-to-end byte-identity vs the Stage-1 fixtures verified by the codegen regen.
- **Not a cap:** a debug `assert` guards nothing away — the overflow silently and correctly routes to
  the heap. `kInlineSyms = 24` is a performance boundary (fast path), never a correctness ceiling.

Optional follow-on: `kMonoAtomInline` 8→4 with the narrowed allocator takes `MonoAtoms` 32→24 and the
term 72→64 (the last of D1's −56%), for atom-light flows.

## Expected cumulative outcome
Stage 1 alone: phase B ~−20–30%, term −12.5%, likely under the 10 GB cap. + Stage 2: term −44%, the
single contraction ~4.5 GB, phase A ~−7–9%. + Stage 3: pathological multiply peak sub-GB. + Stage 4:
the 3ⁿ root cause removed and `dropAqbq47` retired. Each number is filled in from its benchmark as the
step lands — no estimate ships unmeasured.
