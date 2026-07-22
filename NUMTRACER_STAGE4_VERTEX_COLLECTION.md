# Stage 4 — eager quark-gluon vertex collection (retire dropAqbq47)

Goal: keep an internal multi-structure quark-gluon vertex `Σ_a c_a T_a` (T1/T4/T7 of AqbqDirect) as ONE
eager collected object instead of distributing it into 3ⁿ diagrams — exactly, so `dropAqbq47` (which
zeroes T4/T7) is no longer needed. System-wide: helps every flow with internal quark-gluon vertices
(ZAAqbq1/2, ZA3_147, ZA4_147, …), gated behind the existing `DressingCollection -> True` flag so all
other flows stay byte-identical.

## Why it is more than the existing collection

The shipped collection (`$ntDressCollect`, `DSlotOpt`, `dress_collect`, `numeric_value_dressed_netval`;
tested in `tests/test_dpoly.cpp` A–F) collects a propagator numerator `Mq·δ + Z·γ·p̸` — every option's
Dirac structure is **internally contracted** (δ, or a slash `γ·p̸`), NO open Lorentz index. The
vertex structures each carry the **open gluon leg μ**:

- **T1** `γ^μ` — open index on the γ.
- **T4** ~ `p^μ · (δ or slash)` — open index on a *Lorentz vector*, spinor part is δ/slash.
- **T7** ~ `σ^{μν} p_ν` — open index on one leg of a σ, other leg slashed (`dcomm_fs` already exists).

So a vertex slot option is a **(Dirac fragment, open-μ Lorentz factor)** pair. `dress_collect` builds a
concrete Dirac chain per structure-combination and contracts it with `numeric_value_netval`, which
already contracts open Lorentz legs against the net — so the machinery is reused; what is new is (a) a
slot option that can splice an open-μ γ / σ token into the chain AND contribute an open-μ Lorentz
vector into the net, and (b) the front-end decomposition that produces those options.

## Design

### C++ engine (`numeric_contract.hpp`, `dpoly.hpp`)
Extend `DSlotOpt` from `{coeff, dress, slash, vlc}` to also carry an OPEN structure:
```
struct DSlotOpt {
  Cx coeff; std::vector<int> dress;
  // internally-contracted structure (as today):
  bool slash=false; std::vector<pair<double,int>> vlc;
  // NEW — open-gluon-leg structure (μ = the shared vertex axis, a net-known Lorentz id):
  enum Open { None, GammaMu, VecMu, SigmaMu } open = None;
  int openMu = -1;                       // Lorentz id of the gluon leg
  std::vector<pair<double,int>> openVlc; // for VecMu (p^μ) and SigmaMu (σ^{μν}p_ν): the momentum
  // spinor side for VecMu (δ vs slash) reuses the slash/vlc fields above
};
```
`dress_collect` (already enumerates the Cartesian product): per option, in addition to the current
δ/slash splice, emit the open token — `GammaMu → dgamma(openMu)`, `SigmaMu → dcomm_fs(openMu, openVlc)`,
`VecMu → (δ or slash) with an ntVec(openVlc, openMu) factor into the Lorentz net`. The VecMu Lorentz
factor is the one genuinely new routing: it must reach `numeric_value_netval`'s `lor`/`comp`, so either
(i) fold it into the concrete chain as a slashed unit that the net closes, or (ii) pass a per-combination
Lorentz-factor list alongside the concrete chain. Prefer (i) where the structure allows (keeps the
contract closure signature unchanged).

**Checkpoint 1 (foundational, self-contained):** implement `GammaMu` + `SigmaMu` options (open index on
a Dirac token — no new Lorentz routing, reuses `dgamma`/`dcomm_fs`), and add a `test_dpoly.cpp` case:
a slot with options {δ, γ^μ-open, σ^{μν}-open} inside a chain closed by a net carrying μ, collected via
`numeric_value_dressed`, asserted bit-equal to the explicit distributed sum (exactly as tests C–F do).
This proves the open-leg collection mechanism before ANY Wolfram work. `VecMu` (Checkpoint 2) adds the
Lorentz-factor routing.

### Wolfram front-end (`DSL.m`)
Mirror `dressedNumTerm`/`dressedNumDecompose`/`rewriteDressedNums` for the vertex:
- `diracVertexSumQ[p_Plus]` — like `diracNumeratorSumQ` but require the SAME single open Lorentz leg μ
  across all terms (the gluon axis), in addition to the same 2 open spinor indices.
- `vertexNumTerm[t]` — decompose one term into `{din, dout, scalar, openSpec, otherTensors}` where
  `openSpec ∈ {{"gammaMu",μ}, {"vecMu",μ,vlc,spinorSpec}, {"sigmaMu",μ,vlc}}`.
- `vertexNumDecompose[p]` → `(common) * ntVertexSlot[opts, din, dout, μ]`.
- `collectibleDiracSumQ` gains the `|| (diracVertexSumQ && vertexNumDecompose =!= $Failed)` branch.
- `expandVertexSlot[ntVertexSlot[...]]` — inverse (for `redistDiagram`'s cross-check).

### Codegen (`Codegen.m`)
`compileDirac` emits `ntVertexSlot` as a `dtvertexslot` token carrying its options (coeff, dress, open
kind, μ id, vlc) — parallel to how `ntDressedNum` becomes the current slot tokens.

## Verification strategy (exactness is mandatory)
1. **Checkpoint 1 unit test** — collect == distribute, bit-equal, in `test_dpoly.cpp` (no Wolfram).
2. **Redistribution gate** — `redistDiagram` already re-expands a collected diagram and re-analyses with
   collection OFF; extend it to `expandVertexSlot`. The small-D gate in `NumTrace` then cross-checks
   collection == distribution automatically.
3. **Flow-level** — regenerate ZAqbq1_147 (or a cut ZAAqbq) with `DressingCollection -> True` and grade
   the collected kernel against the current distributed kernel via `compare_*` (scale-relative error),
   NOT byte-identity (the collected kernel is a different-but-equal computation). Confirm the net/trace
   count drops toward the dropAqbq47 ≈24-topology scale.
4. **dropAqbq47 retirement** — only after 3 passes: the full-basis collected result must match a
   distributed full-basis reference (small enough flow to generate both) to grading tolerance.

## Order of work
CP1 (C++ open-γ/σ slot + unit test) → CP2 (C++ VecMu Lorentz routing + unit test) → CP3 (Wolfram
`diracVertexSumQ`/`vertexNumDecompose`/`ntVertexSlot` + `expandVertexSlot`) → CP4 (Codegen token) →
CP5 (flow regen + grade + retire dropAqbq47). Each CP is independently testable; CP1–2 need no Wolfram.
