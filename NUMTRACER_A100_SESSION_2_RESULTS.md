# A100 session 2 — results

Answers to the agenda in `NUMTRACER_A100_SESSION_2.md`. One 4-hour interactive allocation on
`n2gpu1211` (A100-SXM4-40GB, CUDA 12.8, 64 cores, 503 GB), 2026-08-11 17:00–21:00.

---

## 0. What had to be rebuilt first, and one deviation from the handoff

The cluster tree was refreshed shortly before the session and **`QCD_Nf2/no_mesons/flows/` was gone**
— it is untracked at `qcd-codes` HEAD (`46b0885`), and there is no Mathematica on this cluster, so it
cannot be regenerated here. HEAD's `no_mesons` driver also no longer matches any available flow set
(66 compile errors: `ZA_integrator::map` argument-list mismatch), and HEAD's YangMills builds but
aborts the flow at the first step (`[Flow aborted early at t = 0.0250]`).

Everything below therefore runs in a **detached worktree at `dd46312`**
(`/pc2/users/c/csdcan02/scratch/Franz/qcd-s1`) — session 1's commit — with:

- `no_mesons/flows/` restored from `numtracer/gpubench/kernels/nf2/flows`, which `vendor.sh`
  vendored verbatim from `no_mesons` **after** the WP1a/WP2b regen. 11 flows, 29,179 generated
  lines. All `numtrace.json` re-pinned `"generated": 1` (ZA4's was 0 — the documented trap).
- `tune_m2A = false` (nf2) and `tune_STI = false` (YM), as session 1.

**Gates pass.** nf2: m2A = 1.1731490534726328 (4.3e-15 vs the handoff's 1.1731490534726379),
Zc(0) = 0.6888524203451546 (1.5e-16), 350 steps. That is the ~1 ulp shift the handoff predicts for
the frame fix, three orders inside the 1e-8 tolerance.

**Deviation, stated plainly:** `dd46312`'s YangMills is the **Logarithmic** variant (WP1b), not the
FocusedLog one session 1 measured — the coordinate type is baked into the generated flow signatures
(`model.hh:63`), so it cannot be switched without regenerating. YM's gate here is therefore
m2A = 6037.240528758882, Zc(0) = 0.9316141347304615, 240 steps — the WP1b values, bit-stable across
every variant below. Angular order is 8 either way, which is the only property YM is being used as a
control for, so the tile conclusions carry; absolute YM wall times are **not** comparable to session
1's table.

---

## 1. The headline, before any knob is turned

Same flows, same config, same GPU. Only DiFfRG changed (`d8e10c83a` → `+ the two fixes`):

| | session 1, LB=128 | session 2, LB=128 | |
|---|---|---|---|
| QCD_Nf2 wall | 22.52 s | **8.15 s** | **2.76x** |

Session 1: `exp5_nf2.tsv` rounds 2–4, 22.51 / 22.54 / 22.50. Session 2: the same binary appears in
four independent A/B sets and measures 8.12 / 8.15 / 8.15 / 8.16 in three of them and 8.31 in the
frame A/B, where it alternated with a 13.9 s variant and ran hotter. Taking the 8.31 instead gives
2.71x; the spread is the honest uncertainty on the last digit, not on the conclusion.

That is the whole of items 3 and 4 arriving at once, and it is larger than every knob session 1
swept, combined.

---

## 2. Item 2 — `ncu --set full`. The measurement no static work substitutes

Nsight Compute, `--set full`, 2 launches per kernel, `--launch-skip 60`, on the rank-5 MDRange
`map` kernels. Reports: `a100_logs/s2/ncu_*.ncu-rep`; summary `a100_logs/s2/ncu_summary.txt`.

### 2.1 Are these kernels fp64-pipe-bound?

**No — but they are not latency-bound either, and the reason has changed.**

Fixed build, `sm__pipe_fp64_cycles_active` as % of peak:

| kernel | fp64 %busy | compute SOL % | memory SOL % | issue slots busy % |
|---|---|---|---|---|
| nf2 ZA4 | 35.8 | 35.7 | 6.9 | 13.4 |
| nf2 ZAAqbq1 | 52.8 | 52.0 | 7.9 | 27.7 |
| nf2 ZAAqbq2 | 52.7 | 51.9 | 7.0 | 25.3 |
| YM ZA4 | 56.6 | 56.0 | 6.7 | 27.4 |

The handoff's laptop roofline estimate — "nf2's big three near ~50% of fp64 peak, YM ZA4 near ~90%"
— is **half right**: the nf2 pair land on ~53%, and ZA4 lower at 36%, but YM ZA4 is 57%, not 90%.

The stated decision rule was "if fp64 utilisation comes back >70%, items 5–6 are pointless." It came
back at 36–57%, so **item 5 was not killed on that test**. It was run anyway, and answered — §6.

### 2.2 Where the frame fix's speedup actually comes from

Fixed vs the same build with the three `device::apply` sites reverted to `tuple_cat`
(divisible tile held ON in both, so this isolates the frame):

| kernel | duration µs (rev → fixed) | local **st** sectors (rev → fixed) | fp64 %busy | long_scoreboard stall |
|---|---|---|---|---|
| nf2 ZA4 | 5.98 → 4.70 (**1.27x**) | 73,267,200 → 20,238,336 (3.6x) | 28.0 → 35.8 | 2.43 → 0.30 |
| nf2 ZAAqbq1 | 3.64 → 1.43 (**2.55x**) | 62,042,112 → 2,709,504 (22.9x) | 18.6 → 52.8 | 3.68 → 0.07 |
| nf2 ZAAqbq2 | 4.65 → 2.16 (**2.15x**) | 68,235,264 → 8,239,104 (8.3x) | 22.8 → 52.7 | 2.74 → 0.07 |
| YM ZA4 | 5.70 → 2.45 (**2.33x**) | 70,385,664 → 3,670,016 (19.2x) | 23.7 → 56.6 | 2.82 → 0.05 |

Three things this settles:

1. **The mechanism is exactly the claimed one.** Achieved occupancy is 12.4–12.5% in *both* builds
   and registers stay at 252–255 — the frame fix buys nothing in occupancy. It removes local-memory
   traffic, and the long-scoreboard stall collapses by 10–50x with it.
2. **YM ZA4's 2.33x on sm_80 reproduces the 4070's 2.323x to three digits.** The frame cost does
   transfer across architectures, as predicted, because you pay it whenever you touch it.
3. **nf2 ZA4 is the outlier at 1.27x, and its local traffic explains why.** It still moves 20.2M load
   and 20.2M store sectors *after* the fix — that is its genuine register spill (1304 B here, §2.3),
   which the frame fix cannot touch. ZAAqbq1 drops to 2.7M and YM ZA4's local **loads go to exactly
   zero**. ZA4 is the one kernel where a spill lever would still pay.

### 2.3 Static confirmation, and a carry-over check that fails

`cuobjdump -res-usage` on `flows/ZA4/src/CT_map_1.cc.o`, LB=128, sm_80, CUDA **12.8**:

| ZA4 | reverted | fixed |
|---|---|---|
| map kernel | 255 reg / **5104 B** stack | 255 reg / **1304 B** stack |
| constant kernel | 48 reg / 3608 B (62.5% occ) | 36 reg / **0 B** (75% occ) |

The handoff §1 asks for one carry-over check before trusting any static number from the laptop:
*nf2 ZA4 → 255 regs / 1864 B at LB=128*. **Registers reproduce exactly; the spill does not** —
5104 B here against 1864 B there. So session-1/laptop *register* counts carry to CUDA 12.8 and
*spill* figures do not. Treat every previously recorded byte count as toolkit-local.

The 3800 B the fix removes from the map frame is the 13 × 272 B interpolator copy plus alignment,
exactly as the commit message claims. The `constant` kernel's frame goes to **zero**, which is the
one place the fix also buys occupancy (62.5% → 75%) — on a kernel worth 0.3–4% of GPU time.

### 2.4 Tile waste, recomputed on this tree

`tilewaste.py --selftest` passes. On the worktree: **nf2 1.255x** masked at LB=128 (the three rank-5
flows at 1.333x each, everything else 1.000x), **YangMills 1.000x** — every tile dimension divides
its extent at angular order 8. Both match the handoff exactly.

The launch geometry confirms the control actually fires. nf2 ZA4, same binary otherwise:

| | block | grid | threads launched |
|---|---|---|---|
| `DIFFRG_DIVISIBLE_TILE=0` (Kokkos) | (32, 4, 1) | (64, 12, 6) | 589,824 |
| `DIFFRG_DIVISIBLE_TILE=1` | (64, 2, 1) | (32, 18, 6) | 442,368 |

589,824 / 442,368 = **1.3333** — the predicted masked fraction, measured off the real launch.

### 2.5 The stall that is left

With `long_scoreboard` gone, the largest remaining issue stall on the fixed build is
`no_instruction` — 3.40 on nf2 ZA4, 0.78–0.94 on the others. On a kernel emitted as 10,684 lines of
straight-line FMA that is instruction-fetch, not data. It tracks kernel size, which makes it the one
metric that argues for *smaller* emitted code rather than better scheduling.

---

## 3. Item 3 — frame fix end-to-end on sm_80

`gpubench-tools/ab.sh`, 4 alternated rounds, round 1 discarded, gates checked every run.

| | reverted | fixed | ratio |
|---|---|---|---|
| QCD_Nf2 | 13.88 / 13.90 / 13.93 s | 8.36 / 8.31 / 8.25 s | **1.67x** |
| YangMills | 5.16 / 5.13 / 5.21 s | 3.34 / 3.32 / 3.35 s | **1.55x** |

Raw: `a100_logs/s2/exp_frame_nf2.tsv`, `exp_frame_ym.tsv`.

Gates: nf2 fixed 1.1731490534726328 vs reverted 1.1731490534726319 (8e-16 apart — the expected FMA
re-contraction); YM bit-identical between the two, since it has 7 interpolators rather than 13 and
nvcc contracted the same way either way.

So the kernel-level 2.2–2.6x lands as **1.55–1.67x on the wall** — the dilution is the host share
the handoff §5 flags as unprofiled (GPU was 64.5% of the Nf2 wall before this fix; after a 2.3x on
the GPU side it is a good deal less, which raises §5's priority rather than lowering it).

One caveat on the `rhs` column in the TSVs: `ab.sh` sums the cumulative `calls:` counters printed on
each residual line, so it counts snapshots, not evaluations. It is stable within an nf2 variant
(626–631 fixed, 701–703 reverted) but varies run-to-run on YM (480 vs 251–255 for the *same*
binary). Do not read it as a work measure.

---

## 4. Item 4 — the tile fix alone, and LaunchBounds re-measured

### 4.1 `DIFFRG_DIVISIBLE_TILE` 0 vs 1 at LB=128, frame fix on in both

| | tile=0 (Kokkos) | tile=1 | ratio | predicted |
|---|---|---|---|---|
| QCD_Nf2 | 10.56 / 10.59 / 10.58 s | 8.16 / 8.17 / 8.16 s | **1.297x** | ~1.20x |
| YangMills | 3.23 / 3.23 / 3.15 s | 3.24 / 3.19 / 3.17 s | **0.998x** | ~1.0x |

The predicted signature lands exactly: a real gain on the code whose angular order is 6, and
nothing at all on the code whose angular order is 8. Gates are **bit-identical** across tile=0/1 on
both codes, as the handoff predicts for a re-tiling that reorders no arithmetic.

nf2 beats its own prediction (1.30x wall against ~1.20x on GPU time) because the fix does two things
at once, visible in §2.4 and in ncu: it stops launching the 1.333x masked lanes **and** it raises
achieved occupancy from 9.4% to 12.5%. Kernel level, tile=0 → tile=1:

| kernel | duration µs | fp64 %busy | achieved occ % |
|---|---|---|---|
| ZA4 | 8.41 → 4.70 (**1.79x**) | 20.5 → 35.8 | 9.45 → 12.48 |
| ZAAqbq1 | 1.87 → 1.43 (1.31x) | 40.8 → 52.8 | 9.43 → 12.49 |
| ZAAqbq2 | 2.78 → 2.16 (1.29x) | 41.3 → 52.7 | 9.41 → 12.48 |

### 4.2 The LaunchBounds sweep, redone

nf2, `DIFFRG_DIVISIBLE_TILE=1`, rounds 2–4 median, with what ptxas did to ZA4:

| LB | wall s | vs 128 | ZA4 regs | spill B | occ | block |
|---|---|---|---|---|---|---|
| 64 | 8.11 | 1.005x | 255 | 1304 | 12.5% | 64 |
| 96 | 8.07 | 1.010x | 168 | 2208 | 18.8% | 64 |
| **128** (shipped) | 8.15 | 1.000x | 255 | 1304 | 12.5% | 128 |
| **192** | **7.81** | **1.044x** | 168 | 2208 | 18.8% | 128 |

YangMills, LB=96 vs 128: 3.19 vs 3.14–3.21 — flat, as in session 1.

**Session 1's sweep is confirmed void.** Its headline — LB=128 an 8–13% local worst, LB=96/64 much
faster — does not survive the tile fix. The whole spread is now 4.4%, and LB=96's advantage over 128
has collapsed from ~8% to 1.0%. The handoff's withdrawal of the `DIFFRG_LAUNCH_BOUNDS=96` arch-key
recommendation is correct, and nothing here reinstates it.

What is left is a clean 2×2, because ptxas splits these four settings into exactly two register caps
and the tile budget splits them into two block sizes:

|  | 64-thread block | 128-thread block |
|---|---|---|
| **255 reg / 12.5% occ** | 8.11 (LB64) | 8.15 (LB128) |
| **168 reg / 18.8% occ** | 8.07 (LB96) | **7.81 (LB192)** |

Register main effect 2.4%, block-size main effect 1.4%, and a clear interaction — LB=192 is better
than either main effect predicts. So occupancy *is* worth something once the tile is fixed, but
2.4%, not the 21% session 1 attributed to it, and it costs 900 B more spill (1304 → 2208 B) to get.

**Recommendation:** LB=192 is a real 4.4% on nf2 and neutral on YangMills, but §6 shows the gain is
*not* the occupancy it looks like, and changing a shipped default on two codes — one of which is
insensitive — is thin evidence either way. LB=128 is defensible today; the withdrawn 96 is not.

---

## 5. Item 7 — per-flow table on the winner, and a correction to §5 of the handoff

`nsys` over the fixed-m2A solve, LB=192 + both fixes, 714 RHS evaluations. Total GPU kernel time
**5.422 s** over 37,099 launches.

| flow | kind | % GPU | total s | calls | avg µs |
|---|---|---|---|---|---|
| ZA4 | map | 41.8 | 2.265 | 714 | 3171.8 |
| ZAAqbq2 | map | 22.6 | 1.226 | 714 | 1716.5 |
| ZAAqbq1 | map | 14.9 | 0.810 | 714 | 1133.9 |
| ZAqbq4 | map | 4.4 | 0.240 | 714 | 335.8 |
| ZAqbq7 | map | 3.5 | 0.191 | 714 | 267.6 |
| ZAqbq1 | map | 2.4 | 0.132 | 714 | 184.3 |
| ZA4 / ZAAqbq2 / ZAAqbq1 | team | 3.2 | 0.176 | 2142 | ~82 |
| ZA, Zc (eta loop) | map | 2.0 | 0.104 | 9048 | ~11 |

The shape is unchanged from session 1: **ZA4 + ZAAqbq1/2 = 79.3%** of GPU time (session 1: 83.1%),
the ZAqbq family 10.3%, the eta loop still 6.34 iterations per RHS and still ~2% of GPU time.
The fixes made everything faster without moving the ranking, so session 1's targeting advice stands.

### The correction: the "host-side block" is mostly GPU wait

The handoff §5 calls host-side work "probably the largest remaining block", on the basis that GPU is
64.5% of the Nf2 wall with "~7.9 s of genuine host compute" remaining, and points at 427,065
`cudaDeviceSynchronize` calls per solve. On this build:

| CUDA API | total | calls |
|---|---|---|
| `cudaStreamSynchronize` | 5.521 s (93.9%) | 35,548 |
| `cudaLaunchKernel` | 0.180 s | 37,099 |
| `cudaMemcpyAsync` | 0.167 s | 52,416 |
| `cudaDeviceSynchronize` | **0.004 s** | **1,533** |

Three consequences:

1. **The 427,065 `cudaDeviceSynchronize` calls are gone** — 1,533 here, 3.6 ms total. Whatever
   session 1 saw (the `async Maps` commit landed between, and it changed exactly the host-side async
   result copy), the "synchronising accessor inside the eta fixed-point loop" hypothesis is dead.
   Do not spend time on it.
2. **`cudaStreamSynchronize` = 5.521 s ≈ total GPU kernel time 5.422 s.** The host is blocked on the
   GPU for essentially the whole time the GPU is busy — so most of what a wall-clock decomposition
   charges to "host" is GPU wait, and it shrank with the GPU. That is why an 8 s wall is possible at
   all: if the ~7.9 s host residue had been genuine compute, 5.4 s of GPU plus 7.9 s of host could
   not fit in 8 s.
3. **The genuine host residue is now ~2.5 s of an ~8 s wall (~31%)**, bounded as
   wall − stream-sync. Still worth the `perf record -g` the handoff proposes — it needs no cluster
   time — but it is no longer "probably the largest remaining block". The GPU is, at ~68%, and
   inside it ZA4 alone is 41.8%.


---

## 6. Item 5 — `minBlocksPerMultiprocessor`, run and closed

The handoff gates this on "only run this if item 2 says the kernels are latency-bound". Strictly it
does not: `long_scoreboard` is ~0.1 on the fixed build. But achieved occupancy sits exactly at the
theoretical maximum and the theoretical maximum is set by registers, so the question was still open,
and there was time. Implemented as the handoff prescribes — a `DIFFRG_MIN_BLOCKS` macro on the
**`dim>1` MDRange specialization only**, so the `dim==1` `RangePolicy` never reaches
`cuda_deduce_block_size` with a `minB` that could drive its block size to 0. The installed header
was restored afterwards; nothing is left patched.

`minB` does exactly what it is supposed to do to registers:

| build | ZA4 regs | spill B | occupancy | wall s (r2–4) | vs LB=128 |
|---|---|---|---|---|---|
| LB=128 (shipped) | 255 | 1304 | 12.5% | 8.15 | — |
| LB=192 | 168 | 2208 | 18.8% | **7.81** | **1.044x faster** |
| LB=128, minB=3 | 168 | 2304 | 18.8% | 8.19 | 0.995x |
| LB=128, minB=4 | 128 | 2688 | 25.0% | 8.95 | 0.911x (**9.8% slower**) |
| LB=128, minB=5 | 96 | 2856 | 31.2% | 9.44 | 0.863x (**15.8% slower**) |

**Two conclusions, and the second is the interesting one.**

1. **Occupancy past 18.8% is a clear net loss, and it gets worse monotonically.** 25% costs 9.8% and
   31% costs 15.8%. The spill bought to get there (1304 → 2688 → 2856 B) is not repaid. `minB` is
   answered and can be closed: do not ship it, and do not sweep it again.

2. **LB=192's 4.4% is not occupancy.** `minB=3` reproduces LB=192's register count (168), occupancy
   (18.8%), *and* its exact launch geometry — both launch grid (32, 18, 6) × block (64, 2, 1),
   verified under ncu — and it is **not faster**; it lands on LB=128. So whatever LB=192 is buying,
   it is ptxas code generation under `__launch_bounds__(192)` versus `__launch_bounds__(128, 3)`,
   not the occupancy that the register count implies. The one visible correlate is 96 B less spill
   (2208 vs 2304).

   One honest caveat on that 4.4%: LB=192 (like LB=96) produces m2A = 1.1731490534726259 where
   LB=128 and every `minB` build produce ...328. Both pass, but LB=192 is on a marginally different
   numerical trajectory, so its wall time is not a pure like-for-like against the others. `minB=3`,
   which *is* bit-identical to LB=128 in gate value, is the cleaner comparison — and it says the
   register cap alone is worth nothing.

This retires the last occupancy hypothesis on the list. After the tile and frame fixes, **occupancy
is not a lever on these kernels in either direction**: raising it via `minB` loses, and the one
setting that appears to win by raising it turns out not to win because of it.

---

## 7. Item 1 — `lambda1L3D` at `-Xptxas -O3`

Submitted as a separate SLURM job on the `normal` partition (`gpubench-tools/s2_item1_o3.sbatch`,
job 34003794, 16 cores / 200 GB) precisely so it could not contend with the wall-clock A/Bs on the
GPU node. Phase 1 builds `with_mesons_3D` as shipped, phase 2 removes the `-O1` pin.

**Phase 1 (shipped `-Xptxas -O1` pin), at -j8: 30 min 08 s wall, peak RSS 42.9 GB.**

That peak is already worth recording: the pin exists because a 30 GB laptop was OOM-killed, and the
comment attributes ~22 GB to *one* `-O3` ptxas. At `-O1` the tree still peaks at **42.9 GB**
(`/usr/bin/time -v` reports the maximum over the process tree, so this is the largest single
process, not a sum). The 30 GB laptop was therefore under-provisioned for this app at *any* ptxas
setting; `-O1` was necessary but not, on its own, an explanation.

---

## 8. Where this leaves the programme

**Settled, do not re-open:**

- The two DiFfRG fixes are worth **2.71x** together on QCD_Nf2 end-to-end on sm_80 (22.50 s → 8.31 s
  at the shipped LB=128), split as ~1.67x frame and ~1.30x tile. Both mechanisms are confirmed at
  kernel level by ncu, not inferred.
- The frame fix transfers across architectures: YM ZA4 2.33x on sm_80 against 2.323x on sm_89.
- The tile fix is worth 1.30x where the angular order is 6 and **exactly nothing** where it is 8.
  `tilewaste.py` predicts which case a code is in, statically and without a GPU. Run it on
  `with_mesons`, `with_mesons_3D` and `finite_T` — the handoff's table says they are all in the
  1.30x class.
- **Occupancy is not a lever here, in either direction** (§6). `minB` closed.
- Session 1's LaunchBounds sweep is void and its `DIFFRG_LAUNCH_BOUNDS=96` recommendation stays
  withdrawn. The residual LB spread is 4.4% and is not occupancy.
- The `cudaDeviceSynchronize`-in-the-eta-loop hypothesis is dead: 1,533 calls, 3.6 ms (§5).

**Still open, in the order I would take them:**

1. **ZA4's remaining spill.** ZA4 is 41.8% of GPU time, sits at 36% fp64 (the lowest of the big
   three), and is the only kernel still moving 20M+ local sectors in each direction after the frame
   fix. Its 1304 B is genuine register pressure. This is now the single largest identified
   inefficiency on the GPU side.
2. **Instruction fetch.** `no_instruction` is the dominant remaining stall on ZA4 (3.40 vs 0.78–0.94
   elsewhere) and it tracks emitted kernel size. That argues for IBP / emission reduction — which is
   also (1)'s lever — rather than for any scheduling knob.
3. **`-Xptxas -O3` on `lambda1L3D`** — §7, finishing as this session ends; the runtime half of the
   comparison still needs a GPU slot.
4. **Host profiling** (`perf record -g`, laptop, no cluster time). Demoted but not dismissed: ~31%
   of the wall, ~2.5 s, and now bounded rather than guessed.
5. **`fastInterpLookups` × launch bounds** (item 6). Untouched. It needs `7ab9cc5fa` since the flag
   no longer exists at HEAD, and its premise — that the flag's cost was unhidden memory traffic that
   might flip at higher occupancy — is weakened by §6 finding no occupancy lever at all. Lowest
   priority, as the handoff had it.

**Two housekeeping items for whoever picks this up:**

- `QCD_Nf2/no_mesons/flows/` is untracked and cannot be regenerated on this cluster. It survives only
  as the `numtracer/gpubench/kernels/nf2` vendored copy. Either commit it or accept that the vendored
  copy is load-bearing.
- HEAD's `no_mesons` driver and HEAD's YangMills do not currently run (§0). Everything here is
  `dd46312`. Reconciling HEAD needs a Mathematica regen, which is laptop work.

**Reproduction.** Builds live under `/pc2/users/c/csdcan02/scratch/Franz/qcd-s1/vacuum/*/build_*`
(`fixed`, `revframe`, `notile`, `lb64`, `lb96`, `lb192`, `mb3`, `mb4`, `mb5`). Raw TSVs, ncu reports,
nsys report, res-usage dumps and both patched header variants are in `a100_logs/s2/`. Scripts added:
`gpubench-tools/s2_build_variants.sh`, `s2_build_minb.sh`, `s2_run_item4.sh`, `s2_item1_o3.sbatch`,
`ncu_summary.py`. The installed DiFfRG headers were restored to the shipped state and verified
byte-identical to the `DiFfRG_KT` working tree — nothing is left patched.
