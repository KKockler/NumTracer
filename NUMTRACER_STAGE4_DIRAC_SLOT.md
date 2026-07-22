# Stage 4 (redesigned) — the general collected Dirac slot: ANY number of open legs

**Supersedes `NUMTRACER_STAGE4_VERTEX_COLLECTION.md`.** That doc (and the committed CP1/CP2 C++)
modelled a vertex slot as a single open leg carried by one token (`γ^μ` / `σ^{μν}` / `p^μ`, the `Open`
enum). Ground truth from a real `AqbqDirect147` vertex (probe against `za3_147`, 2026-07-22) shows the
model is too narrow on **two** independent axes, and the AAqbq basis breaks it on a third:

1. **Options are multi-token Dirac chains, not single tokens.** The open gluon leg is always a `γ^ρ`
   closed by the gluon propagator/projector; the *variety* is in the internally-contracted chain
   **before** it. Measured: T1 = `γ^ρ`; T4 = `p̸₁·γ^ρ` (slash **then** open γ); T7 =
   `[γγ·v − γγ·v]·γ^ρ` (a 2-term commutator prefix, then the open γ).
2. **`VecMu` (open index on a vector) does not occur** in `AqbqDirect147` — every leg attaches through
   a γ. It is a real *possibility* (below), just not the common case CP2 assumed.
3. **A vertex can have MANY open legs.** The full `{A,A,qb,q}` (AAqbq) vertex has **two** open gluon
   legs per vertex; a general n-gluon quark vertex has **n**. Some AAqbq structures put both legs on a
   **non-Dirac** tensor, e.g. `g^{μν}·δ_spinor` (both legs on a metric, trivial spinor) — the k=2
   analogue of CP2's `p^μ`. The single-`int` `openMu` cannot express this at all.

The optimal target is therefore: **arbitrary open-leg count k ≥ 0, legs on either the Dirac side or the
Lorentz side.** This turns out to be *natural*, not a stretch — see feasibility below.

## The general object

A **collected Dirac slot** is a runtime-coefficient-weighted sum of Dirac structures that all share:
- the same spinor-in `din` and spinor-out `dout` (one quark-line position), and
- the same **set** of open Lorentz legs `{μ₁,…,μ_k}` (k ≥ 0), each closed by the surrounding net.

| k | physical object | today |
|---|---|---|
| 0 | propagator numerator `Mq·δ + Z·p̸` | `ntDressedNum` (shipped, `DSL.m`) |
| 1 | single-gluon quark vertex `Aqbq` (T1/T4/T7) | distributed → 3ⁿ (Stage 4 target) |
| 2 | two-gluon quark vertex `AAqbq` | distributed (worse) |
| n | general | — |

This is **one mechanism**. k=0 is the propagator collection already shipped; k≥1 is the vertex. The
redesign unifies them under one head (`ntDiracSlot`) and one C++ type (`DSlotOpt`), so there is exactly
one code path to validate.

Each **option** (summand) = `numeric coeff × dressing atoms × structure`, where **structure** =
- a **Dirac-token chain** `din→dout` (γ / slash / γ5 / σ / commutator …) whose *free-Lorentz tokens*
  are open legs on the Dirac side, **times**
- a set of **Lorentz-net factors** (vector `p^μ`, metric `g^{μν}`, …) carrying the remaining open legs.

The union of the two token/factor sets' free Lorentz ids **must equal** the slot's shared leg set
`{μ₁,…,μ_k}` for every option (the net closes a *fixed* set of legs regardless of structure choice).

## Feasibility — why "any k" is almost free

The C++ contraction engine **already** closes an arbitrary number of free Dirac legs against the net:
`numeric_dirac` folds a closed spinor loop into a `Factor` over **all** its free gluon legs (guarded at
f ≤ 15), and `contract_factors`' greedy variable-elimination contracts that tensor against the net's
matching legs. This is exactly how the existing 3-gluon-vertex **quark triangle** — many free gluon legs
per loop — is contracted today. Nothing in the engine assumes one open leg. The "single leg" assumption
lives only in three places, all replaceable:
- CP1/CP2 `DSlotOpt` (`openMu` is one `int`, `open` is one enum) → replace with a token list + net-factor list.
- The plan's Wolfram predicate ("the SAME single open Lorentz leg μ") → a leg **set**.
- Codegen (`openMu` single) → the toks' free legs + the net factors.

So the general design is the token-chain design (already agreed) plus: (i) net-factor list instead of the
single `VecMu` vector, (ii) leg **set** instead of one μ throughout the front-end. No new contraction math.

## C++ engine — the general `DSlotOpt`

Replace the `Open` enum + `slash`/`vlc`/`openMu`/`openVlc` fields with:

```cpp
struct DSlotOpt {
  Cx coeff{1,0};
  std::vector<int> dress;               // dressing-atom ids (product = the runtime coefficient)
  std::vector<network::DFac> toks;      // Dirac-token chain spliced in place of the slot (din→dout);
                                        // its free-Lorentz tokens are this option's Dirac-side open legs
  std::vector<network::Elem> netFacs;   // extra Lorentz-net factors for this option (open legs on
                                        // vectors/metrics/…); appended to every net term before contraction
};
```

`dress_collect` per chosen option: splice `opt.toks` into the concrete `DiracNet`, and stash `opt.netFacs`
to append to the Lorentz net (via the closure, exactly as CP2 routes `VecMu` vectors — generalised from
"a vector" to "any Lorentz Elems"). `numeric_value_netval` then closes every free leg — Dirac-side and
net-side — automatically. The `tr(1)=4` collapse restoration (`nCollapsed`) is unchanged.

This **subsumes** every prior case:
- propagator δ → `toks={}`, `netFacs={}`; slash → `toks={dslash(vlc)}`.
- CP1 `GammaMu` → `toks={dgamma(μ)}`; `SigmaMu` → `toks={dcomm_fs(μ,vlc)}`.
- CP2 `VecMu` → `netFacs={Vector(μ,p)}` (+ optional `toks={dslash}` for the spinor side).
- real T4 → `toks={dslash(p1), dgamma(μ)}`; T7 → the two commutator chains as two options.
- AAqbq `g^{μν}δ` → `toks={}`, `netFacs={Metric(μ,ν)}`; two open γ's → `toks={dgamma(μ), dgamma(ν)}`.

## Wolfram front-end (`DSL.m`) — general k

Unify with the shipped propagator collection rather than parallel it:
- `collectibleDiracSlotQ[p_Plus]` — every term is a Dirac structure with the **same** `{din,dout}` and
  the **same set** of open Lorentz legs (k≥0). k=0 recovers `diracNumeratorSumQ`; k≥1 is the vertex.
- `diracSlotTerm[t, legs]` — decompose one term (after pulling common non-Dirac factors) into
  `{din, dout, scalar, structureSpec}` where `structureSpec` is a **general Dirac-chain + net-factor**
  serialisation: an ordered token list `{{"gamma",μ}|{"slash",vlc}|{"sigma",legA,legB}|{"g5"}…}` plus a
  net-factor list `{{"vec",μ,vlc}|{"met",μ,ν}…}`, whose combined free legs equal `legs`.
- `diracSlotDecompose[p]` → `(common) * ntDiracSlot[opts, din, dout, legs]`.
- `expandDiracSlot[ntDiracSlot[...]]` — inverse (for `redistDiagram`'s small-D cross-check).
- Register `ntDiracSlot` in `tensorQ`/`labelsOf` (**returns `Join[legs,{din,dout}]`** — the legs are real
  Lorentz labels the net contracts), `spinorLabelsHead` (`{din,dout}`), and wherever `ntDressedNum` is
  listed. Eventually `ntDressedNum` becomes the k=0 alias (or is folded away).

To obtain each option's canonical chain, pull common non-Dirac factors out, then `Expand` **only** the
Dirac part into monomial chains (the vertex sum is small; this is the same expand `splitColourGroupsInv`
already does, but the result is *kept collected* as slot options instead of distributed into diagrams).

## Codegen (`Codegen.m`)

`compileDirac` treats `ntDiracSlot` like `ntDressedNum` (a `dtslot`); `diracSlotStrBody` emits each
option as `DSlotOpt{coeff, {dress}, {toks…}, {netFacs…}}` — the toks as `dgamma/dslash/dcomm_*/dg5`
builders (free-leg ids = `ids[μ]`), the net factors as `network::Elem{…}` literals. Multi-leg falls out
of serialising the lists. Extend the same head-lists as `ntDressedNum` (diracIn/out, orderDiracFacs
token count, the leak guards, `splitColourGroupsInv`).

## New staging

Each checkpoint independently testable; exactness mandatory (grade on numbers where the emitted set
changes, byte-identity where it must not).

- **R0 — general `DSlotOpt` (token chain + net factors), C++ only.** Replace the `Open` enum with
  `{toks, netFacs}`; rewrite `dress_collect` to splice `toks` and append `netFacs`. Migrate tests C–H to
  the new form. **New unit test:** a k=1 option that is a *multi-token* chain matching the real T4/T7
  shape (`{dslash(p1), dgamma(μ)}`), collected == distributed, bit-equal. Retires CP1/CP2's enum.
  *(No Wolfram; L1 `ctest` green.)*
- **R1 — multi-leg engine validation, C++ only.** A **k=2** slot (two open legs) — one option with two
  open γ's (`{dgamma(μ), dgamma(ν)}`), one with a metric on the net (`netFacs={Metric(μ,ν)}`, δ spinor),
  one with a mixed vec+γ — closed by a 2-leg net, collected == distributed over random points. Proves
  `numeric_value_netval` closes an arbitrary leg count through the slot. *(No Wolfram.)*
- **R2 — Wolfram decomposer, general k.** `collectibleDiracSlotQ` / `diracSlotTerm` / `diracSlotDecompose`
  / `ntDiracSlot` / `expandDiracSlot`; the general chain+net serialiser; head-list registration.
  **Round-trip test (Wolfram):** `expandDiracSlot@diracSlotDecompose[v] ≡ v` (modulo leg relabelling)
  on a real k=1 (`AqbqDirect147`) vertex sum AND a real k=2 (`AAqbq`) vertex sum, captured pre-distribution.
  *(This is the milestone this session targets.)*
- **R3 — unify the propagator path + gate.** Route `collectibleDiracSumQ`/`distributeQ`/`rewriteDressedNums`
  through `ntDiracSlot` (k=0 = today's `ntDressedNum`); keep `DressingCollection` gating. Regen the
  dressed kernels; grade byte/numeric-identical (k=0 must not move any value).
- **R4 — codegen token, general.** `compileDirac` emits `ntDiracSlot` with multi-leg `toks`/`netFacs`.
  L2 on the small flow; byte-identity where the fast path is untouched.
- **R5 — flow grade + retire the `AqbqDirect1` restriction.** Regen `za3_147` (k=1) and a cut `AAqbq`
  (k=2) with the FULL basis via the collected slot; grade collected-vs-distributed (scale-relative);
  confirm the net/trace count drops toward the collected scale; retire the internal-vertex
  classical-only approximation.

**This session:** R0 → R1 → R2 (reach the Wolfram round-trip), per the agreed scope. R3–R5 follow.

### Progress

**R0 + R1 — general `DSlotOpt` (token chain + net factors) — LANDED 2026-07-22.**
`DSlotOpt` is now `{Cx coeff, vector<int> dress, vector<network::DFac> toks, vector<network::Elem> netFacs}`
— the `Open` enum and `slash`/`vlc`/`openMu`/`openVlc` fields are gone. `dress_collect` splices each
option's `toks` into the concrete chain and appends its `netFacs` to every net term; the closure's
byte-identical fast path is `netFacs`-empty. This retires the committed CP1/CP2 single-token model and
subsumes it. All `test_dpoly` cases C–H migrated to the new form (still exact). New cases:
- **I (R0):** multi-token k=1 options `{γ^μ}`, `{p̸₁,γ^μ}`, `{p̸₁,γ^μ,p̸₂}` (open leg first / last / middle),
  collected == distributed, **0.00e+00** over 5000 pts × two parities (T1 and multi-token T4 each nonzero).
- **J (R1):** a **k=2** slot (two open legs {μ,ν}) with options `{γ^μγ^ν}` (both legs Dirac-side),
  `g^{μν}·δ` (`netFacs={Metric(μ,ν)}`, both legs net-side), and `p̸·γ^μ · p^ν` (mixed) — closed by a
  2-leg net, all three nonzero (|tr| 43/43/16), collected == distributed to ≤7e-15. Proves
  `numeric_value_netval` closes an arbitrary open-leg count through the slot.

All 35 non-codegen tests green.

**R2 — Wolfram decomposer, general k — LANDED 2026-07-22.**
`DSL.m` gains `ntDiracSlot[opts, din, dout, legs]` (registered in `tensorQ`/`scalarQ`/`labelsOf`
[= `Join[legs,{din,dout}]`]/`spinorLabelsHead`), and:
- `openLorentzOf[t]` / `colourLabelsOf[e]` — a term's open Lorentz legs = free indices that are neither
  spinor nor colour (the gluon axes).
- `diracSlotSumQ[p]` — a Plus whose EXPANDED terms each carry Dirac structure, the same 2 open spinor
  indices, and the same NON-EMPTY open-leg set (k≥1; k=0 stays on the `ntDressedNum` path).
- `diracSlotDecompose[p]` → `(commonColour)(commonScalar) · ntDiracSlot[opts, din, dout, legs]`, each
  option `{residualDressing, structureProduct}` with the structure kept WHOLE (Dirac chain × its
  Lorentz-net factors, e.g. the gluon projector on the open leg); the toks/netFacs split is deferred to
  codegen (R4). `Expand` first so an inner Dirac Plus (a σ commutator) splits into monomial options.
- `expandDiracSlot` — inverse (for `redistDiagram`).

Design choice: the slot's open legs = the **sum's free Lorentz legs** (the external gluon axis reached
through the common propagator), NOT the internal γ index — so no per-term relabelling is needed and the
common projector factors out cleanly. Validated by round-trip (`expand@decompose ≡ Expand[sum]`):
- synthetic k=1 (open leg via a projector) and k=2 (Dirac legs / a metric / mixed vec+γ) — exact
  (`tests/gen/test_diracslot_roundtrip.wls`).
- the **three real `za3_147` `AqbqDirect147` vertex sums** captured pre-distribution — exact, legs
  `{v_i}`, 6/14/9 options each (`Expand` distributes the T7 commutator into monomial options; each slot
  is ONE collected trace vs 3ⁿ diagrams). `tests/gen/probe_vertex_sum.wls`.

**R3 — collection gate — LANDED 2026-07-22.** `collectibleDiracSumQ` now also accepts a vertex sum
(`diracSlotSumQ && diracSlotDecompose =!= $Failed`), so `expandBridges`/`distributeQ` keep it eager and
`rewriteDressedNums` rewrites it to `ntDiracSlot` (k=0 propagator numerators still take the `ntDressedNum`
path); `redistDiagram` expands `ntDiracSlot` too. Verified: `collectibleDiracSumQ=True` and
`rewrite→ntDiracSlot` for the synthetic k=1 and k=2 sums.

**R4 — codegen token, general k — LANDED 2026-07-22.** `compileDirac` treats `ntDiracSlot` like a
`dtslot`; `diracSlotStrBody` splits each option's whole structure into `DSlotOpt{coeff, {dress}, {toks},
{netFacs}}` — the Dirac chain ordered via a new `orderOpenChain` (an OPEN din→dout walk) and emitted as
`dgamma`/`dslash`/`dcomm`/`dg5` tokens, the Lorentz factors emitted as `network::Elem` literals via a new
`builderInvElem` (projectors carry their env `Base`/`Inv` atom, so the open leg closes exactly as the
distributed diagram's). Internal legs (the γ↔projector bridge) keep their own ids and momenta because
`allLabels`/`momentumOf` recurse into the slot via `Cases[Infinity]` — no fresh-id allocation. The
`ntDressedNum` propagator emitter was migrated to the same `{toks, netFacs}` form (its `netFacs` always
empty). Two bugs found & fixed en route: emit bare `Elem` (the gen TU has `using namespace
numtracer::network`, no `network::` alias); and `din/dout` must follow the chain ORIENTATION (a `Sort`
reversed the token order → wrong trace).

**R5 — full-basis flow grade — LANDED 2026-07-22 (k=1).** Regenerated `za3_147` (the FULL `AqbqDirect147`
1/4/7 quark-gluon vertex) with collection ON: the vertex now folds to ONE `ntDiracSlot` per leg instead
of 3ⁿ diagrams. `compare_za3_147_num` grades the collected kernel: **with ZAqbq4=ZAqbq7=0 it reduces to
the FORM-validated struct-1 kernel at 6.3e-11 rel (pure round-off, threshold 1e-9)**, and differs with the
real 4/7 dressings — i.e. the collected full-basis vertex is physics-exact. This single regen exercised
BOTH the new `ntDiracSlot` vertex path AND the reformatted `ntDressedNum` propagator emitter (the quark
propagators are k=0 slots). All 35 non-codegen tests stay green.

**The measured diagram win** (`tests/gen/probe_za3_147_count.wls`, analysis only): `za3_147` with the full
147 vertex is **6 diagrams collected vs 221 distributed — 36.8× fewer** (net top terms 6 vs 437). The
collected kernel is byte-different but physics-identical, at a 37× smaller diagram set.

**BUT generation TIME gets WORSE, not better** (`tests/gen/probe_za3_147_timing.wls`) — an important,
counterintuitive finding:

| phase | collection ON (6 diag) | collection OFF (221 diag) |
|---|---|---|
| NumTrace (analysis) | 0.08 s | 3.21 s |
| generator compile (clang++) | 2.7 s | 3.1 s |
| generator **run** (reduce+rebase+lower) | **26.9 s** | **0.34 s** |
| TOTAL generation | **31.6 s** | 13.3 s |
| emitted kernel size | 178 KB | 239 KB |

The 36.8× fewer diagrams makes the Mathematica *analysis* 39× faster and the kernel ~25% smaller, but the
**generator RUN is 79× SLOWER (26.9 s vs 0.34 s)**, so total generation is **2.4× slower with collection**.
Cause: a collected diagram is ONE big `DPoly` (all structure×dressing combinations as polynomial
coefficients); reducing + CSE/Horner-lowering that monolith is costly, whereas the distributed path emits
~221 small independent traces that the shipped **cross-trace CSE + sub-term dedup** collapse to a compact
shared program in 0.34 s. Collection front-loads the 3ⁿ enumeration into one `DPoly` the deduper cannot
factor across — it trades the two runtime-lowering optimisations for the diagram-count/RAM/kernel-size win.

**Implication.** Stage 4's payoff is diagram count, Mathematica-analysis time, emitted **kernel size**, and
(expected) contraction **RAM** — NOT generator wall-time. On a small compile-light flow like `za3_147` the
DPoly-lowering cost dominates and collection is a net time LOSS.

### Root cause of the generator-run regression (profiled, `NT_GEN_PROFILE=1`)

```
[num] phase A: 6 distinct traces, 0 cached, table 0.0 MB, 0.0 s (W=32)   ← contraction is FREE
[num] phase B+lower: 6 nets in 4 groups, window 6, 26.9 s (W=32)          ← ALL the cost is here
```
The 26.9 s is **entirely phase-B lowering**, and it is **serial**:
- In `fold_groups_streaming` (`trace_fold.hpp`) the per-net **fold** is parallel (`parallel_flat`, W=32),
  but the **`sink` — the `to_genprog(acc)` lowering — runs on the calling thread in the drain loop**,
  once per group in ascending order. That serial order is deliberate: every group lowers into ONE shared
  `GlobalEnv g`, whose interning is what makes the kernel small (178 KB) and reproducible.
- With collection there are only ~6 nets / 4 groups, each a **big DPoly** (the summed structure×dressing
  combinations), so lowering is ~4–6 serial `to_genprog` calls on 32 cores — ≤6-way, 26 idle cores.
- Compounding it: the collected DPoly is a bigger lowering job than the 221 distributed traces, because
  distribution keeps them small + **dedup'd + cross-net-interned** (`best_into`/`GlobalEnv` sharing folds
  the 221 to a compact program in 0.34 s), whereas the collected sum is one large unique polynomial the
  CSE cannot factor across.

**The tension is fundamental:** the small kernel comes from SHARING (serial interning into one
`GlobalEnv`, optionally `FusedStream` cross-CSE), and sharing forces serial lowering. The two things the
user liked — small kernel — and disliked — slow generation — are two sides of the same coin.

### CORRECTION (measured 2026-07-22): the bottleneck is the CONTRACTION, not the sink

Implementing the parallel sink (optional mutex on `GlobalEnv::intern` + parallel per-group `to_genprog`)
and regenerating `za3_147` moved phase B **not at all** (26.9 s → 28.1 s). `NT_GEN_POLYSTATS=1` shows why:
the collected DPolys are **small** (monos 6 / 18 / 113 / 1568; SSA ≤ 4717) — lowering is *fast*. The 28 s
is the dressed **contraction** (`dress_collect` enumerating the quark-triangle option combinations), and
it runs **inside phase B's fold**, not phase A: because no trace recurs (`nReused=0 → nCache=0`), the
traces are contracted **on demand** in `fold_net`. The fold is parallel over nets, but there are only ~6
nets and ONE dominates (its `dress_collect` does ~37 `numeric_value_netval` contractions **serially**).
So the earlier "serial sink" reading was wrong: the sink is cheap; the serial cost is `dress_collect`.
The parallel-sink change was reverted (correct but it optimises the wrong phase; shipping an unvalidated
optimisation violates the measure-everything discipline).

**The real fix: parallelise `dress_collect`'s combination enumeration** (each combination is an
independent `contract` → accumulate; a flat parallel work list with a serial merge). The wrinkle is
nesting — `dress_collect` runs *inside* the phase-B fold's `parallel_flat`, so either (a) make the
dressed fold serial and let `dress_collect` use the cores (few dressed nets, so ≈ full parallelism, no
nesting), or (b) expose each collected net's combinations to phase A's flat work list (the "matrix-of-
DPoly" contraction, plan 4b) so the existing phase-A parallelism covers them. `numeric_contract.hpp`
can't include `trace_fold.hpp` (circular; the net-builder TUs forbid its threads), so this needs a small
shared parallel primitive or the fold restructured. This is the correct next optimisation.

### (b) — execution-ready implementation plan (the real fix: restore cross-net dedup)

Root cause recap (measured): the collected net's `dress_collect` contracts EVERY structure×dressing
combination separately, serially per net (za3_147: 26.9 s, one dominant net), whereas the distributed
path contracts only the DISTINCT traces (5–7× cross-net dedup) over phase A's flat parallel list
(0.34 s). (b) routes the combinations through that same table so both the dedup and the parallelism
re-apply. `dress_collect` (the C++ enumerate-and-contract) leaves the hot path entirely.

The single structural change: **a sub-term's scalar must carry a dressing monomial**, so a collected
combination becomes an ordinary trace (its concrete Dirac chain, deduped) times a `(Cx × ∏dress atoms)`
scalar, and the net's DPoly is assembled by the fold. Steps, with the exact touch points:

1. **Codegen — enumerate at emit time (`compileDirac` / `splitColourGroupsInv`, `Codegen.m`).** For an
   `ntDiracSlot` net, expand the Cartesian product of slot options into one **plain** sub-term per
   combination: `{concreteDiracNet, ln, dressMonoAtoms, Cx}` — the concrete chain is the fixed tokens +
   the chosen options' `toks`, the net gets the options' `netFacs` appended, and the scalar is the
   product of the options' `coeff` × dressing-atom ids. Emit `ntDressedCore` no more; these are ordinary
   `{sdn, sln}` sub-terms that dedup in the trace table exactly like the distributed ones.
2. **Sub-term scalar carries the dressing monomial (`Codegen.m` sub-term build ~L2148 + the emitted
   `dsc`).** Today `dscv : vector<vector<Cx>>`. Add a parallel `vector<vector<vector<int>>>` of dress-atom
   id lists (one per sub-term), or fold the atom ids into the scalar literal. Emitted alongside `dsc`.
3. **Fold path — MPoly trace → DPoly by a dressing scale (`trace_fold.hpp` / the `scale` lambda in the
   phase-B emission).** The dressed backend's `T` becomes **`MPoly`** (plain deduped traces), and the
   `scale(d, MPoly&&)` lambda returns a **DPoly** = `fromMPoly(trace)` scaled by `Cx` and keyed by the
   sub-term's dress monomial (`out.add(dmono, trace*Cx)`). The accumulator stays `DPoly`; only the trace
   type and the scale signature change. `fold_net`/`fold_groups_streaming` are already templated on
   `<P>` for the accumulator — the change is that trace-type ≠ accumulator-type for this path, so add a
   second template param (trace `TP`, accumulator `P`) or a thin dressed-fold wrapper.
4. **Retire the C++ enumeration from the hot path.** `numeric_value_dressed_netval` / `dress_collect`
   stay for the unit tests and as the small-D reference, but the generator no longer calls them.
5. **Validation.** (i) `test_dpoly`: a new case asserting the assembled DPoly (plain traces × dress
   monomials) bit-equals `numeric_value_dressed`. (ii) Regen `za3_147`, confirm `[cse]` now reports the
   deduped distinct-trace count (≈ the distributed scale) and phase B drops to **~0.3–1 s**; grade with
   `compare_za3_147_num` (must still reduce to struct-1). (iii) Confirm the emitted kernel stays the
   compact DPoly form (parametric in dressings — the 178 KB win is preserved because the DPoly is still
   assembled per net, only its contraction is deduped).

Expected: `za3_147` generation **31.6 s → ~4–5 s** (below distributed's 13.3 s), bounded phase-B RAM (W
live traces), kernel unchanged-small. Serves both the time and the ZAAqbq RAM north star. This is a
focused contraction-pipeline change (~2 files of real logic + the fold generalisation) that must land as
one validated unit — do NOT half-apply it (a sub-term scalar half-carrying a dress monomial silently
drops dressings).

#### SIMPLER realization (codegen-only — the one to build)

The scalar/fold change above is avoidable. Instead of ONE dressed sub-term per net carrying the whole
multi-option slot (so C++ `dress_collect` enumerates the 3ⁿ combinations serially at contraction time),
**emit each combination as its OWN dressed sub-term whose slot has a SINGLE option** (the chosen
combination: its concrete chain in `sdch`, that combination's coeff+dressing in `sdsl`). Then:
- the sub-term key `(sdn, sln, sdch, sdsl)` dedups combinations **across nets** in the existing trace
  table (the measured 5–7×) — and within a net the ~37 combinations are distinct anyway;
- phase A contracts the distinct combinations over its **flat parallel** list (each is now a trivial
  1-combination `numeric_value_dressed_netval` = one contraction, so no serial 3ⁿ loop anywhere);
- the net's DPoly is assembled by the **existing** phase-B fold (each combination's trace is a
  one-dressing-term DPoly, summed) — no scalar-type or fold-type change.

So (b) reduces to a **codegen-only** change: where `compileDirac` emits `ntDressedCore[chain, slots]`
(one sub-term), instead expand the slot's Cartesian product and emit **N single-option dressed
sub-terms**. The C++ `dress_collect` stays (it now runs with one option per call — the deduped,
parallel, on-demand path handles the 3ⁿ). This is the tractable implementation; validate as in step 5.

#### LANDED 2026-07-22 — codegen-only expansion (the run regression is FIXED; compile is the new bottleneck)

Implemented (Wolfram-only, `Codegen.m`), no C++ change:
- `dressedSlotStrBody` / `diracSlotStrBody` now return the **bare option-string list** (not the wrapped
  `DSlot{…}`); `compileDirac`'s dressed return carries the slots **structured** (`ntDressedCore[chain,
  slots]`, `slots` = list-over-slots of option-lists).
- `emitNumericGenerator` expands each dressed branch's Cartesian product (`Tuples[slots]`) into **one
  single-option sub-term per combination** — same `DiracNet{}`/rest/scalar/chain, a `DSlot` list with
  ONE `DSlotOpt` per slot. Non-dressed branches unchanged (byte-identical).
- **The missing half the plan didn't call out:** the singletons only reach phase A if they are *cached*.
  Phase A contracts `[0, nCache)` and Codegen.m defaulted `nCache = nReused` — but every expanded
  combination is a **distinct singleton** (`nReused ≈ 0`), so the reused-only default left ALL of them
  to phase B's fold, which is parallel over NETS not traces ⇒ the one dominant net serialised thousands
  of contractions (unchanged 27 s). Fix: **default `nCache = NSUB` for dressed flows** (each collected
  combination is individually tiny — see table). `NT_GEN_MEMO_MAX` still dials it back for RAM.

Built in THREE codegen-only pieces (miss any one and it regresses), measured on `za3_147` (full
`AqbqDirect147` vertex), collection ON:

| phase | OLD collected | expand + `nCache=NSUB` | **+ pool-interning (final)** | distributed (OFF) |
|---|---|---|---|---|
| NumTrace | 0.08 s | 0.08 s | 0.09 s | 3.21 s |
| generator source | — | 9.5 MB | **1.06 MB** | — |
| generator **compile** | 2.7 s | 13.5 s (39 units) | **3.7 s (8 units)** | 3.1 s |
| generator **run** | **26.9 s** | 2.1 s | **2.5 s** | 0.34 s |
| — phase A | 6 traces, 0 cached | 12101 cached (8.1 MB), 1.9 s | 12101 cached (8.1 MB), 2.3 s | — |
| — phase B+lower | 26.9 s | 0.0 s | **0.0 s** | — |
| TOTAL generation | ~31.6 s | ~19.3 s | **~8.9 s** | 13.3 s |
| emitted kernel | 178 KB | 178 KB | **178 KB** | 239 KB |

The three pieces:
1. **Expand** each combination into a single-option sub-term (above). Moves the 3ⁿ enumeration off C++
   `dress_collect` onto the trace table — but only if the singletons reach phase A (piece 2).
2. **`nCache=NSUB` default for `hasDressed`** — the combinations are all distinct singletons
   (`nReused≈0`), so the reused-only default left them to phase B's net-skewed serial fold (run stayed
   27 s). NSUB puts them on phase A's flat W=32 list ⇒ run 26.9 → 2.1 s; phase B+lower 0.0 s.
3. **Pool the chain + option columns.** Expansion makes both columns ~99.9% redundant — measured on
   `za3_147`: **84672 `DSlotOpt` emissions but only 35 DISTINCT**, and 12101 chains but 5 distinct. So
   emit the pools (`chp`, `optp`) ONCE + per-sub-term INDEX arrays (`sdchR`, `sdslR`), and rebuild
   sdch/sdsl in main O(nSub) — the exact hash-consing the `sidx`/`dsc` tables already use. Source 9.5 MB
   → 1.06 MB ⇒ compile 13.5 → 3.7 s. (Without piece 3 the expansion trades the run regression for a
   COMPILE regression; with it, neither.)

**Correctness:** value-exact — each single-option sub-term is exactly one term of the old `dress_collect`
odometer sum (same coeff × dressing monomial, same `nCollapsed` tr(1)=4). `flow_za3_147_num` reduces to
struct-1 at round-off and differs with real 4/7 (PASS). Kernel renumbered (value-identical). The **k=0
propagator flow `za3_num` regenerates BYTE-IDENTICAL** (`flow_za3_num` PASS) both with and without the
pooling — the strongest signal the common path is untouched. Pooling is a pure source-level change (the
emitted kernel is identical to the pre-pooling expanded kernel).

**Result: collection is now the FASTEST path on `za3_147` — ~8.9 s total gen vs distributed 13.3 s vs
old-collected 31.6 s — while keeping the 178 KB kernel.** The generator-RUN regression is gone
(26.9 → 2.5 s) AND the compile stayed flat (2.7 → 3.7 s). Stage 4's payoff is no longer *just* diagram
count / kernel size / RAM — on this flow it now wins wall-time too.

**ZAAqbq measured (dense proxy, `gen_zaaqbq1_small_numeric.wls` — the only suite flow with real sub-term
dedup).** Regenerated with the full change: **byte-identical kernel, `flow_zaaqbq1_small` PASS.**
- generator source **1.01 MB**, compile 3.3 s, run 1.5 s, total gen **8.65 s**;
- phase A: **1180 distinct traces, all cached, table 11.9 MB**, peak RSS **~125 MB** (trimmed to 42 MB);
- redundancy **8.3× PRESERVED** (dedup-OFF ref = 9792 traces → 1180 deduped), so the expansion did NOT
  dissolve the cross-net reuse here (unlike za3_147, where every combination was unique per net).

So the feared **phase-A RAM blow-up from `nCache=NSUB` is a non-issue** on this flow — the collected
combinations are individually tiny (11.9 MB for 1180) and, where redundancy is high, `nReused ≈ NSUB`
anyway so NSUB caches only marginally more. The pool-interning keeps the source ~1 MB. `NT_GEN_MEMO_MAX`
remains the dial if a bigger flow's phase-A table ever grows.

**Caveat — the ULTIMATE target is still untested:** the small proxy uses `ZAqbq1` (struct-1 internal
vertex), so its vertex slots have few options. The genuine payoff — the FULL-basis `AqbqDirect147`
collected ZAAqbq (retire the `AqbqDirect1` classical-only restriction) — would generate FAR more
combinations per net (za3_147-scale × a 4-point vertex), and it is NOT yet wired. That switch is where
phase-A RAM under `nCache=NSUB` must be re-measured. But the mechanism (expand + NSUB + pool) is now
proven correct and cheap on both a high-combination flow (za3_147, 12101, no reuse) and a
high-redundancy dense flow (ZAAqbq1 small, 8.3× reuse).

**Remaining:** wire `AqbqDirect1 → AqbqDirect147` for ZAAqbq/ZA and grade the full-basis k=2 flow (+
re-measure phase-A RAM there); regen the other committed dressed fixtures for repo consistency (all
become renumbered but value-exact — mechanism validated on k=0 + k=1 + the dense proxy, not yet swept);
update `tools/redundancy.py`, which parses inline `dch`/`dsl` per sub-term and no longer sees the pooled
index arrays (the generator's own phase-A / dedup-OFF lines report the redundancy directly meanwhile).

### Fix options (superseded by the correction above — kept for the design record)
1. **Parallelise the sink into thread-local `GlobalEnv`s, then merge** — lower the 4 groups' DPolys
   concurrently, remap+dedup their programs at the end. Up to ~min(nGroups, cores)× on this phase; keeps a
   (near-)shared final kernel. The real work is a correct, cheap env-merge (remap instruction refs, dedup
   constants). Best long-term answer; a non-trivial C++ change in `trace_fold.hpp` + `gen.hpp`.
2. **Parallelise CSE, serialise only interning (RECOMMENDED — the split IS clean).** Confirmed by reading
   `lower_into`: the heavy `gdetail::best_into` (CSE/Horner) writes ONLY the thread-local `RBuilder w`; the
   sole shared-`g` mutation is the id interning (`g.var_id/inv_id/dr_id`) while *building* the monomials,
   which is cheap hash inserts. So: put a mutex on `GlobalEnv::intern`, pre-size `progs`, and run the
   per-group sink (`progs[gi]=to_genprog(acc,g,ro)`) via `parallel_flat` — the CSE runs fully concurrent,
   only the interning serialises. Expected ≈ min(nGroups, cores)× on phase B (za3_147: ~4 groups → the
   26.9 s → ~7 s, making collection *faster* than distribution's 13.3 s). **Caveats:** (a) interning under
   a lock is order-nondeterministic → the env is RENUMBERED vs a serial run (value-identical, not
   byte-identical) — fine since collected kernels are value-graded, but it also perturbs the
   NON-collected flows' byte fixtures, so it must be validated by the `compare_*` value grades across all
   35 flows, and either gated (only when nGroups·polysize is large) or accepted as a global renumber; (b)
   it slightly serialises every flow's interning (mutex overhead) — measure it stays negligible.
3. **Independent per-group lowering (no shared env)** — trivially parallel, but a BIGGER kernel (loses
   cross-group interning). Directly trades the kernel-size win for generation speed; measure both.
4. **Gate collection to RAM-bound flows** — collect only where distribution's RAM forces low `nB`
   (`ZAAqbq`); distribute the compile-light flows (`za3_147`). Sidesteps the regression but forfeits the
   kernel-size win on the non-collected flows.

**Recommendation.** Measure the RAM/`nB`/time picture on a real `ZAAqbq` regen first (the north-star flow —
its distributed lowering is itself the bottleneck there, so collection may WIN on time too, not just RAM).
Then pursue (2)→(1) to reclaim the phase-B parallelism while keeping the small kernel. This reframes W1 as
a **RAM/size** lever whose generator-run cost is a real, now-understood regression with a clear fix path.

**Remaining.** A k=2 flow-grade (a flow with an INTERNAL two-gluon `AAqbq` vertex — the mechanism is
already validated by engine test J + the synthetic k=2 decompose/gate, only an end-to-end flow is
untested); and switching the `ZAAqbq`/`ZA` setups from `AqbqDirect1` to `AqbqDirect147` via collection
(the north-star RAM/speed payoff), with a performance check that collection keeps the net count near the
collected scale.

## Ground-truth evidence (probe `tests/gen/probe_vertex_sum.wls`, captured pre-`expandBridges`)

`AqbqDirect147` vertex, open gluon leg `v1` via `ntTransProj[p1,v1,ρ]` (the common gluon propagator):
- T1 (`ZAqbq1`): `γ^ρ` — 1 token, open ρ.
- T4 (`ZAqbq4`): `p̸₁ · γ^ρ` — slash then open γ (2 tokens).
- T7 (`ZAqbq7`): `[γ^{μ}γ^{ν}·v·v − …] · γ^ρ` — a 2-term commutator prefix, then open γ.
The distribution happens in `expandBridges` (via `distributeQ`, `DSL.m`) — extending
`collectibleDirac*Q` keeps the sum eager; `rewriteDressedNums` then rewrites it to the slot token;
`splitColourGroupsInv` (codegen) force-expands γ-bearing **Pluses** but not single tokens, so the slot
survives to `numeric_value_dressed_netval`.
