# NumTracer — performance summary

Generated numeric kernels vs their FormTracer (FORM) reference, on this machine
(`-O3 -march=native`, `bench_aqbq147` over 200k random points). `num/FORM` is the ratio of
kernel evaluation time; correctness is the max relative error against FORM.

## Current headline (2026-07-19): the numeric backend is FASTER than FORM

`bench_aqbq147`, 200k points, `taskset`-pinned:

| struct  | FORM ns | numeric ns | num/FORM | was  | rel-err vs FORM |
|---------|---------|------------|----------|------|-----------------|
| ZAqbq1  | 1237.9  | **1184.0** | **0.96x** | 2.90x | 3.04e-09 |
| ZAqbq4  | 1130.0  | **1117.6** | **0.99x** | 3.02x | 9.74e-11 |
| ZAqbq7  | 2040.0  | **1262.8** | **0.62x** | 1.84x | 2.78e-09 |

Three changes landed to get there, in order of size:

1. **Multi-term denominator cancellation** (`divThroughPolyAtoms`, `numeric/mpoly.hpp`) — the big one.
   `divThroughMonomialAtoms` only ever cancelled an inverse atom whose denominator was a single
   *monomial*, so a shifted line `k = l − q` (with `k²` multi-term) kept its `1/k²` atom even when the
   numerator carried a factor of that same `k²` — which Dirac-trace numerators routinely do. Terms are
   grouped by their atom multiset and each multi-term denominator is trial-divided in; an exact
   division (vanishing remainder, same relative tolerance as the noise prune) drops the atom. On
   ZAqbq1_147 Mq-in: 384 of 1,168 trial divisions are exact, monomials 13,269 → 4,832, SSA
   33,775 → 7,649, runtime 3021 → 1102 ns. **Accuracy improves** (one division and a pile of
   cancelling terms are gone). Generation gets *faster*. `NT_GEN_NO_POLYDIV=1` disables.
2. **`CrossTraceCSE` / `trace_all`** (`codegen/gen.hpp`, `numeric/numeric_contract.hpp`) — all trace
   polynomials lower through one shared CSE program. 2.64x → 2.47x on ZAqbq1_147 Mq-in. Opt-in
   (`"CrossTraceCSE" -> True`): it costs ~2.2x on the consumer's `-O3` compile.
3. **Coefficient snapping** (`snap_coeff`, `numeric/numeric_contract.hpp`) — 2.90x → 2.64x, free.

This supersedes the 1.35/1.37/1.05 row below, which had been unchanged since the initial commit and
predates `tests/gen/gen_qcd_aqbq147_numeric.wls` — it never measured this flow.

### Emitted-size effect per flow (`kernels.hh`, before → after the denominator cancellation)

The win is concentrated where shifted-line propagators meet matching numerators — quark and ghost
loops. Pure-gauge and small unit-test flows are untouched.

| flow | before | after | |
|---|---|---|---|
| ZAqbq1_147 / ZAqbq4_147 / ZAqbq7_147 | 1.92M / 1.81M / 2.01M | 431K / 389K / 451K | **0.22x** |
| ZAqbq1, ZAcbc | 74,902 / 50,931 | 20,482 / 13,944 | **0.27x** |
| ZAAqbq1Small (+Ref) | 88,213 / 89,341 | 42,660 / 43,000 | 0.48x |
| Zq | 4,375 | 2,693 | 0.62x |
| ZA4, YM_ZA4 | 179,714 / 163,486 | 137,266 / 125,572 | 0.76x |
| ZA4_147 | 5,122,566 | 4,261,663 | 0.83x |
| ZA3, ZA3_147 | 21,249 / 232,498 | 18,501 / 211,930 | 0.87x / 0.91x |
| ZA, EtaPiL, Discdirac, Glu_quark, Pion_quark, unit flows | — | unchanged | 1.00x |

`ZA3_147` runtime is unchanged at 1.01x num/FORM (`bench_147`) — a pure-gauge 3-point has little to
cancel. `ctest` is 36/36 including `codegen_regen` after regenerating every flow above.

## Symbolic dressing collection (default-on, no gate)

Keeping a dressed propagator-numerator sum eager as one collected `DPoly` trace (instead of
distributing into `2^D` diagrams) is now a **pure win** — faster than the previous GlobalCollect
path on both runtime and generation, so it stays default with no gate.

**(historical — see the current headline above; these ratios were measured on a flow definition that
predates `tests/gen/gen_qcd_aqbq147_numeric.wls` and do not reproduce today.)**

Fundamental-colour quark–gluon vertex `ZAqbq{1,4,7}_147` (struct 1 / 4 / 7):

| metric                | GlobalCollect (old) | collection, unfixed | **collection (current)** |
|-----------------------|---------------------|---------------------|--------------------------|
| runtime, num/FORM     | 1.42 / 1.45 / 1.10  | 1.74 / 1.72 / 1.35  | **1.35 / 1.37 / 1.05**   |
| generation (struct 1) | 27 s                | 45 s                | **27 s**                 |
| rel-err vs FORM       | —                   | —                   | 2.5e-7 / 5.1e-9 / 3.5e-8 |

What made it a win (both fixes in `mathematica/Codegen.m`, source-only — non-dressed kernels
regenerate byte-identical):

- **Runtime:** the imaginary-vanishing probe (`numericImagProbeRealQ`) is now dressed-aware. It
  previously failed to compile against a dressed `fill()` (missing the `dr_<id>` dressing-atom
  args), so every dressed kernel was conservatively kept complex and computed-then-discarded its
  imaginary half. It now stubs and passes the dressing atoms, confirms `Im≈0`, and emits the
  lossless **RePart** double kernel (re/im split; the dead imaginary halves are DCE'd by the
  compiler) — exactly what the distributed path already got.
- **Generation:** the dressed generator's main TU compiles at `-O1` (env `NT_GEN_MAIN_OPT`) — it
  runs once (~1 s) so `-O2`'s extra cost on the `DPoly`/`dch`/`dsl` templates was wasted (28.7 s
  vs 17.5 s, identical run). The `lc/dc` net-CSE also now runs for the dressed path.

`to_genprog(DPoly)` prunes round-off **per dressing channel** (each kinematic poly against its own
max), not against one global tolerance — a channel reweighted small at generation time can dominate
at runtime once a regulator suppresses the others, so a global tolerance could silently delete its
genuine terms.

## Other flows

All committed numeric kernels match their FORM / equivalence oracle — the full ctest
suite passes, including the per-flavour / per-component / disconnected flows
(`flow_flavour_split`, `flow_flavour_ingroup`, `flow_gluon_condensate`, `discdirac_num`,
`flow_sigl_dirac`). For indicative ratios of the larger gluon kernels (quark box `ZA4` ≈ 1.15×,
three-gluon `ZA3` 1/4/7 within a small factor of FORM), see the per-flow `compare_*` drivers and
`tests/refshim/`. The generation pipeline (the disposable build-time generator program) is
described in the repo's top-level `README.md` (the *Architecture* section).

## Reproducing

```bash
# regenerate a flow (needs wolframscript + FunKit, and FORM for the YM/QCD tensor bases):
wolfram -script tests/gen/gen_qcd_aqbq147_numeric.wls
# runtime + rel-err vs FORM (not a ctest — compile by hand):
g++ -O3 -march=native -std=c++20 -Iinclude -Itests/gen -Itests/refshim \
    tests/refshim/bench_aqbq147.cpp -o /tmp/b && /tmp/b
# the full numeric gate:
ctest --test-dir build --output-on-failure
```
