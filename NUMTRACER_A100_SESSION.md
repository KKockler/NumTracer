# A100 measurement session (WP5) — self-contained instructions

**Audience:** a fresh Claude Code session (or Franz) running **on the Noctua2 cluster**
(PC2 Paderborn, A100 nodes). This document assumes no prior context — read it fully before
running anything.

**Mission:** runtime-measure, on a real A100 (sm_80), the four knobs that could not be
resolved on the RTX 4070 Laptop (sm_89, fp64 1:64 — its slow fp64 pipe hides register/spill
effects that matter on a 1:2 datacenter part). This is a **measurement session**: change only
the experiment knobs described here, record numbers, commit nothing. Defaults get decided
afterwards, together, back on the laptop.

The four experiments, in recommended execution order:

| # | Experiment | Knob | Cost |
|---|---|---|---|
| 1 | Baseline + physics gate + `fastInterpLookups` A/B | runtime config flag | 2 builds, 4+ solves |
| 2 | Per-flow GPU profile | `nsys` on baseline binary | 2 profiled solves |
| 3 | Launch-bounds sweep {96, 128, 192} | `-DDIFFRG_LAUNCH_BOUNDS=N` compile flag | +2 builds per app |
| 4 | `__noinline__` size-gate sweep {off, 300, 500, 800} | `noinline.py` rewrite of `kernels.hh` | +3 variant trees per app |

Then: a winner-combination run, and the IBP go/no-go numbers (§9).

---

## 1. Context: where things stand (2026-08-10)

Two production fRG codes, both DiFfRG-based, GPU-accelerated via Kokkos/CUDA:

- **YangMills** (`qcd-codes/vacuum/YangMills/SP`, binary `build/YangMills`): 5 flows
  (ZA4, ZA3, ZAcbc, ZA, Zc). FocusedLog coordinates. The solve is **86.9% ZA4**.
- **QCD_Nf2** (`qcd-codes/vacuum/QCD_Nf2/no_mesons`, binary `build/QCD`): 11 flows.
  GPU time (4070 nsys, fixed-m2A solve, 714 RHS): **ZA4 39.1%, ZAAqbq2 26.9%,
  ZAAqbq1 17.1%** (= 83% in three kernels), ZAqbq family 9.3%, rest small.

Kernel composition (SASS attribution, 4070): nf2 ZA4 = **77% NumTracer-emitted arithmetic**
(giant polynomial traces); ZAAqbq1/2 = **57–58% interpolator** (coordinate-transform
`log`/`log1p`/`asinh` per lookup); regulators only 4–11% everywhere.

GPU round 2 (already landed, all physics-gated): production flow regen (~1.08x), position
cache in `QuadratureIntegrator::map()` (neutral on 4070, expected A100 upside), k-only lookup
hoisting (−450 fp64/kernel, +1.5–2.5%), and **`fastInterpLookups`** — a tabulated inverse
coordinate transform (runtime config flag, default false) worth **1.13x e2e on YM** on the
4070. Cumulative on the 4070: ~1.25x YM, ~1.10x Nf2. Full numbers:
`NumTracer_bckp/NUMTRACER_GPU_ROUND2.md` and
`NumTracer_bckp/numtracer/gpubench/results/baseline_2026-08-09/NOTES.md`.

**Why the A100 changes the picture.** fp64 is 1:2 on A100 vs 1:64 on Ada → ~40x more fp64
throughput on a fixed problem. Arithmetic stops hiding everything else: registers, spill,
occupancy, and memory behaviour become visible. Static cross-compilation for sm_80/sm_90
(done on the laptop, `ptxas`/`cuobjdump` cross-target fine) predicts:

- **Occupancy is poor**: nf2 ZA4/ZAqbq1/ZAqbq4 sit at **12%** on sm_80 (LaunchBounds<128>,
  ~168–255 regs).
- **The `__noinline__` size gate**, currently OFF in all production kernels: static sweep on
  nf2 ZA4 showed threshold **500** (the sm_89-tuned default) leaves 255 regs / 12% occupancy
  on sm_90, while **300** gives 0 spill / 168 regs / 18% occupancy at +2.6% fp64 ops, and
  **200** gives 25% occupancy at +5% ops. sm_80 spill numbers were not captured — that is
  what this session measures. On the 4070 the gate measured −4.8% wall (real but small);
  on a 1:2 part the occupancy/spill side of the trade is far heavier.
- `-maxrregcount` is **INERT** (launch bounds take precedence — verified byte-identical
  output for all values). The LaunchBounds *value* is the only register knob;
  LB ≥ 512 collapses to 64 regs with 3–5x worse spill — do not test it.

## 2. What must be on the cluster (run these rsyncs FROM THE LAPTOP first)

Three trees. The flow state in `qcd-codes` is **not fully in git** (regenerated
`flows/**` are working-tree state) → rsync, don't clone. SSH alias `noctua2`
(or `n2login1`) is in `~/.ssh/config`. Pick a target dir on the parallel FS
(`$PC2PFS/<project>/...`), call it `$T` below:

```bash
T=noctua2:/scratch/<project>/<user>/a100_session   # adjust!

# 1. DiFfRG (branch feature/fast_interpolators, must include commits
#    20c5a5f4a "fast interpolator lookups" and cd588376d):
rsync -av --exclude build --exclude bin --exclude _deps --exclude Testing \
      /mnt/data/Documents/Uni/Code/DiFfRG_KT/ $T/DiFfRG_KT/

# 2. The two production codes (exclude build trees and bulky artifacts):
rsync -av --exclude build --exclude 'build_*' --exclude TBCache --exclude figs \
      --exclude '*.nsys-rep' --exclude '*.sqlite' \
      /mnt/data/Documents/Uni/Code/qcd-codes/vacuum/YangMills/ $T/YangMills/
rsync -av --exclude build --exclude 'build_*' --exclude TBCache --exclude figs \
      --exclude '*.nsys-rep' --exclude '*.sqlite' \
      /mnt/data/Documents/Uni/Code/qcd-codes/vacuum/QCD_Nf2/ $T/QCD_Nf2/

# 3. gpubench tools + reference results (small):
rsync -av /mnt/data/Documents/Uni/Code/NumTracer_bckp/numtracer/gpubench/tools/ \
      $T/gpubench-tools/
rsync -av /mnt/data/Documents/Uni/Code/NumTracer_bckp/numtracer/gpubench/results/baseline_2026-08-09/ \
      $T/ref-results/

# 4. This file:
rsync -av /mnt/data/Documents/Uni/Code/NumTracer_bckp/NUMTRACER_A100_SESSION.md $T/
```

**The NumTracer engine / Wolfram is NOT needed on the cluster.** The noinline sweep does not
regenerate anything: `gpubench-tools/noinline.py` reproduces the generator's size gate as a
pure text rewrite of a vendored `kernels.hh` (same transformation as `eff_decor` in
`gen.hpp`, threshold as a parameter). Verified on the laptop to reproduce the gate's output.

## 3. Cluster environment + builds

Discover, don't assume — but the shape on Noctua2 is (no `-A` account flag needed, the
user's default association covers it; adjust partition via `sinfo` if `gpu` is wrong —
there is prior SLURM history on this cluster, old job scripts in `~/` are a good template):

```bash
# This Claude session should ALREADY be running inside an interactive GPU job — the user
# starts it that way, because a session cannot drive `srun --pty bash` itself (each Bash
# call is a fresh non-interactive shell). Verify at session start:
echo $SLURM_JOB_ID && nvidia-smi -L   # expect a job id + A100-SXM4-40GB
# If instead you find yourself on a LOGIN node inside an `salloc` allocation
# ($SLURM_JOB_ID set but no GPU visible): prefix every build/run/profiling command
# with `srun`, which executes it on the allocated node.
# If neither is set, stop and tell the user to re-spawn you inside:
#   srun -p gpu --gres=gpu:a100:1 --cpus-per-task=32 --mem=64G -t 8:00:00 --pty bash
# (generous -t: when the job dies, this session dies with it)

module avail cuda gcc cmake   # then e.g.:
module load CUDA/12.x GCC CMake
nvidia-smi   # note ECC state and clocks in the results log
```

Rules:
- **Never compile on login nodes**; do everything inside the job. Compute nodes have
  128 cores / 256+ GB — `-j32` is fine there (the laptop's `-j8` rule does not apply).
- Do **not** build DiFfRG's test suite and do not run ctest — waste of time here.

### 3.1 DiFfRG

If a DiFfRG install already exists on the cluster (likely — check `~/.local/share/DiFfRG`,
old `CMAKE_PREFIX_PATH`, previous job scripts): update its **source** to the rsynced
`DiFfRG_KT/DiFfRG` state and rebuild `libDiFfRG` + install. The round-2 changes are headers
**plus `src/DiFfRG/common/configuration_helper.cc`** (reads the `fastInterpLookups` flag), so
a header-only `cmake --install` is NOT enough — the library must recompile.

If no install exists: `DiFfRG_KT/install.sh` is the superbuild (builds deal.II etc. —
hours; do it once with `-j32` on a compute node).

CUDA arch: configure with the GPU visible (Kokkos autodetects) or set
`-DKokkos_ARCH_AMPERE80=ON` explicitly. **Verify after the first app build:**
`cuobjdump -lelf build/QCD | head` must show `sm_80`.

**Trap:** `cmake --install` preserves source mtimes. After reinstalling headers, consumers
may silently recompile **zero** translation units. After any DiFfRG change:
`find <install>/include -name '*.hh' -exec touch {} +`, then check the app build log
actually shows `Building CXX/CUDA object` lines.

### 3.2 Apps

```bash
cd $T/YangMills/SP
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=<DiFfRG install> \
      # + whatever the existing cluster setup used (mirror old job scripts if present)
cmake --build build -j32          # builds app + flow kernels; do NOT build test targets
# same for $T/QCD_Nf2/no_mesons  (binary name: QCD)
```

Do **not** run any `numtrace`/regeneration target — the rsynced `flows/**` are the
validated round-2 state. If CMake tries to regenerate (missing Wolfram will make it fail
or skip), that's fine as long as `flows/*/kernels.hh` are untouched; check
`git status`-style mtimes or just verify the solve gate (§5) passes.

## 4. Measurement protocol (applies to EVERY experiment)

1. **Fixed workload**: time the *gate configs*, not the adaptive production solves:
   - Nf2: `no_mesons/run_gate/parameter.json` (fixed-m2A solve; on the 4070 it does
     exactly **714 RHS evaluations** — a controlled workload).
   - YM: `SP/parameter.toml` (the committed production one — the solve is short, ~17–20 s
     on the 4070). Record the step count from the log; if two variants differ in step
     count, compare **time per step**, not total wall.
2. **Alternate variants**: never A A A B B B. Run A B A B A B (binaries in separate build
   dirs make this cheap), **discard round 1** (warmup: FS cache, CUDA context, JIT),
   compare **ratios per round**, report the median of ≥2 kept rounds. A100s hold clocks
   far better than the laptop, but node-sharing and thermal effects still exist. Log
   `nvidia-smi --query-gpu=clocks.sm,temperature.gpu --format=csv` before each round.
3. **Physics gate before any timing claim** (per binary, once): the solve must reproduce
   - YM (FocusedLog): `m2A = 6037.244194877355`, `Zc(0) = 0.9316318995037505`
   - Nf2 (run_gate): `m2A = 1.1731490534726277`, `Zc(0) = 0.6888524203451526`
   to **≤1e-8 relative** (expect e-11..e-13: different GPU/CPU/compiler shifts last ulps;
   bit-identity with the laptop is NOT expected and NOT required). The final line of
   `output.log` prints both values. Reference H5s (laptop, pre-round-2 HEAD, physics-equal
   to e-11): `$T/ref-results/{ym,nf2}_head_dd46312_output.h5`.
   ⚠ The laptop's `SP/build/parameter.toml` was left in an experimental state (its last
   run **diverged**) — use the committed `SP/parameter.toml`, not a build-dir copy.
4. **Invocation**: apps read `parameter.json`/`parameter.toml` from the CWD, or take
   `-p <file>`. Run each solve from its own run dir (outputs land in CWD):
   `cd run_gate && ../build/QCD` — this is exactly how `run_gate/` was produced.
5. **Static SASS caveat**: both interpolator index paths (analytic + table) are always
   compiled — the `fastInterpLookups` branch is runtime. Static fp64 counts of interp code
   are therefore inflated ~2x; on this session, grade on **runtime**, use
   `cuobjdump -res-usage` only for registers/spill deltas.
6. Record everything in `$T/A100_NOTES.md` as you go (format in §10).

## 5. Experiment 1 — Baseline + `fastInterpLookups` A/B

The flag is `/discretization/fastInterpLookups` (bool, default false), read once at
interpolator construction. `run_gate/parameter.json` has it **false**;
`run_fastinterp/parameter.json` has it **true**; for YM, toggle the line in a copy of
`parameter.toml` (`fastInterpLookups = true` under `[discretization]`) — or override on the
CLI if supported (`./QCD --help` shows the syntax; editing the file is the safe path).

Per code: gate-check both flag states once, then alternate OFF/ON per §4.2.

**Questions this answers:** (a) A100 baseline wall for both codes (the denominator for
everything else); (b) does the 4070's 1.13x YM / 1.02x Nf2 from the table transform grow or
shrink on a 1:2-fp64 part? (Expectation: interp share *drops* when arithmetic gets 40x
faster relative to memory, but the ZAAqbq1/2 kernels are 57–58% interp — Nf2 could gain
*more* than on the laptop. Genuinely open.)

Deliverable: `wall(OFF)/wall(ON)` per code + gate values.

## 6. Experiment 2 — per-flow GPU profile (runs while other builds compile)

On the baseline binary (flag ON — the intended production state):

```bash
cd run_gate
nsys profile -o nf2_a100 --stats=true ../build/QCD
# kernel table: nsys stats --report cuda_gpu_kern_sum nf2_a100.nsys-rep
```

Same for YM. Record the per-flow GPU-time shares table (kernel names contain the flow
names). This is the direct A100 analogue of the 4070 numbers in §1 and the input to the
IBP decision (§9). Also extract: total GPU time vs wall (CPU share of the solve!), and
the number of RHS evaluations.

## 7. Experiment 3 — LaunchBounds sweep {96, 128, 192}

The knob is the `DIFFRG_LAUNCH_BOUNDS` macro (`DiFfRG/common/kokkos.hh`, default 128,
documented override `-DDIFFRG_LAUNCH_BOUNDS=N`). It parameterizes the Kokkos launch policy
templates instantiated in **app** TUs, so rebuilding the app with the define is sufficient
(no DiFfRG rebuild):

```bash
for LB in 96 192; do
  cmake -B build_lb$LB -DCMAKE_BUILD_TYPE=Release <same opts as build> \
        -DCMAKE_CXX_FLAGS="-DDIFFRG_LAUNCH_BOUNDS=$LB" \
        -DCMAKE_CUDA_FLAGS="-DDIFFRG_LAUNCH_BOUNDS=$LB"
  cmake --build build_lb$LB -j32
done
```

(128 = the existing `build`.) Verify the define took: `cuobjdump -res-usage` register
counts must differ between build dirs for the big kernels. Gate each binary once, then
alternate all three per §4.2. Priority: **Nf2 first** (the occupancy-starved code), YM if
time permits.

Decision input: best LB per code; a change from 128 needs **≥2%** e2e to be worth it.

## 8. Experiment 4 — `__noinline__` size-gate sweep {off, 300, 500, 800}

The headline experiment. Currently every production `kernels.hh` is fully inline
(`grep -c noinline` = 0 across all 16 flows). The gate out-of-lines trace functions with
more than N SSA statements as `KOKKOS_FUNCTION __attribute__((noinline))`, which is the
**only** lever that removes the giant-kernel spill on datacenter parts (no LaunchBounds
value does it on sm_90; sm_80 to be measured here).

### Variant creation (seconds, no Wolfram):

```bash
cd $T/QCD_Nf2/no_mesons
python3 $T/gpubench-tools/noinline.py flows/ZA4/kernels.hh --report   # see what gates where

for MIN in 300 500 800; do
  cp -a . ../no_mesons_ni$MIN            # full source copy (flows are #included at compile time)
  for f in flows/*/kernels.hh; do
    python3 $T/gpubench-tools/noinline.py "$f" -o "../no_mesons_ni$MIN/$f" --min $MIN
  done
  find ../no_mesons_ni$MIN/flows -name kernels.hh -exec touch {} +   # cp -a kept mtimes
done
```

Same pattern for `YangMills/SP` → `SP_ni$MIN`. Then configure+build each copy into its own
`build/` (fresh build dirs — full compile anyway). Note: only flows with functions over the
threshold change; small flows compile identically — that's fine.

**Verify before trusting a variant:** `grep -rc noinline flows/*/kernels.hh` — nf2 ZA4 must
show hits at 300/500 (7 functions >500 lines, `tr0` alone ~2,280) and the diff vs the
original must touch **only decorator lines** (`diff` should show
`KOKKOS_INLINE_FUNCTION` → `KOKKOS_FUNCTION __attribute__((noinline))` and nothing else).

### Measurement

Per variant: physics gate once (results were **identical** to the inline build on the
4070 — expect the same), record `cuobjdump -res-usage` regs/spill for the ZA4 kernel, then
alternate {off, 300, 500, 800} per §4.2. Nf2 first, then YM.

Expected from static analysis: 300 best (0 spill, 18%+ occupancy, +2.6% ops), 500
marginal, 800 ≈ off. If 300 wins clearly, a **200** point is worth adding
(25% occupancy at +5% ops — the trade tips toward occupancy on a 1:2 part).

Decision input: threshold minimizing Nf2 e2e; ≥2% to justify changing the shipped default.
Note for later (do NOT do it on the cluster): if a non-500 threshold wins, the default in
`gen.hpp` becomes arch-keyed, and the power-table + caller-side trace-call-hoisting
questions (closed as "premise void" while the gate was dead) reopen.

## 9. Winner combination + IBP go/no-go numbers

1. Build one binary per code with the winning LB **and** winning noinline threshold **and**
   `fastInterpLookups` at its better setting; alternate vs the §5 baseline. This is the
   headline number: **A100 e2e speedup of the session**.
2. Re-run the nsys profile (§6) on the winner. From it + the composition percentages in §1,
   compute the IBP ceiling on A100:
   `gain ≈ Σ_flow share(flow) × NTshare(flow) × 0.3` (0.3 = optimistic fraction of trace
   arithmetic an integration-by-parts reduction could remove; 4070 numbers: nf2 ZA4
   NTshare = 0.77, ZAAqbq1/2 ≈ 0.16–0.22, YM ZA4 ≈ 0.20). On the 4070 this gave ~15% YM /
   ~35% Nf2 ceilings. If the A100 numbers keep the Nf2 ceiling ≥20%, IBP is GO for a
   future round; otherwise NO-GO. **Record the inputs, don't start IBP.**

## 10. Recording + bringing results home

Append to `$T/A100_NOTES.md`, same style as `ref-results/NOTES.md`. Per measurement:
date, node/GPU (`nvidia-smi -L`), clocks/temp, module versions, binary (build dir + flags),
config file, gate values, per-round walls, ratio table. Plus the two nsys kernel tables
and a final summary block:

```
| experiment | variant | Nf2 ratio | YM ratio | regs/spill (nf2 ZA4) |
```

At the end, from the laptop:
`rsync -av $T/A100_NOTES.md $T/*.nsys-rep <laptop>:/mnt/data/Documents/Uni/Code/NumTracer_bckp/numtracer/gpubench/results/a100_2026-08-11/`

## 11. Condensed trap list

- **cmake --install preserves mtimes** → consumers skip rebuild; touch installed headers,
  verify TUs actually compile (§3.1).
- **cp/rsync -a preserve mtimes** → touch copied `kernels.hh` before building variants.
- **Do not edit the extended lambdas** in `quadrature_integrator.hh`: `if constexpr`
  INSIDE an extended device lambda miscompiles under nvcc (crash or silent zero-RHS).
  Nothing in this session requires touching them.
- **Static SASS of interp code is ~2x inflated** (both index paths compiled); runtime only.
- **`-maxrregcount` is inert** (launch bounds win) — don't waste a build on it.
- **LB ≥ 512 collapses** to 64 regs / 3–5x spill — excluded by design.
- Adaptive step counts can differ between variants → compare per-step/per-RHS when counts
  differ (§4.1).
- Don't build test suites; don't compile on login nodes; `-j32` on compute nodes is fine.
- If a solve diverges, the physics gate failed — stop and diagnose the config (wrong
  parameter file is the usual cause, see §4.3), don't time it.
- This is measurement only: no engine/DiFfRG source edits beyond the documented knobs,
  no commits, no default changes.
