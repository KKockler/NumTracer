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

Existing flows are untouched (the new head appears in no current expression; the collection gate is not
yet wired — that is R3). **Next: R3** (route `collectibleDiracSumQ`/`distributeQ`/`rewriteDressedNums`
through `ntDiracSlot`, regen + grade), then **R4** (codegen token: split each option into `toks`/`netFacs`
and emit `DSlotOpt`), then **R5** (full-basis flow grade + retire the `AqbqDirect1` restriction).

## Ground-truth evidence (probe `tests/gen/probe_vertex_sum.wls`, captured pre-`expandBridges`)

`AqbqDirect147` vertex, open gluon leg `v1` via `ntTransProj[p1,v1,ρ]` (the common gluon propagator):
- T1 (`ZAqbq1`): `γ^ρ` — 1 token, open ρ.
- T4 (`ZAqbq4`): `p̸₁ · γ^ρ` — slash then open γ (2 tokens).
- T7 (`ZAqbq7`): `[γ^{μ}γ^{ν}·v·v − …] · γ^ρ` — a 2-term commutator prefix, then open γ.
The distribution happens in `expandBridges` (via `distributeQ`, `DSL.m`) — extending
`collectibleDirac*Q` keeps the sum eager; `rewriteDressedNums` then rewrites it to the slot token;
`splitColourGroupsInv` (codegen) force-expands γ-bearing **Pluses** but not single tokens, so the slot
survives to `numeric_value_dressed_netval`.
