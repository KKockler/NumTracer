# Per-flavour and per-component dressings

This is the most physics-heavy tutorial; the [Glossary](../getting_started/glossary.md) defines
*dressing*, *flavour*, *isospin*, and *condensate*. Its runnable artifacts are the Mathematica
generators (`tests/gen/gen_flavour_split.wls`, `gen_flavour_ingroup.wls`, `gen_gluon_condensate.wls`)
and their committed kernels — validated by the ctests named below — rather than a standalone
`Tutorials/` C++ program.

The flows so far dressed every quark the same way — one propagator dressing $Z_q$ for the
whole colour/flavour multiplet. Physics often needs more: the **u** and **d** quarks dressed
differently (broken isospin), or a colour component dressed on its own (a gluon condensate).
NumTracer supports this two ways, and the choice matters for how big the kernel gets.

| approach | how | when |
|---|---|---|
| **per-flavour split** | one diagram (and one dressing) per flavour | a handful of fields, dressed independently; isospin is *not* a symmetry |
| **group-diagonal dressing** | one trace, a $\delta$ tagged with named dressings on selected components | dress selected components *within* an SU($N$) group without splitting the trace |

The dividing line: folding a closed flavour loop as an SU($N_f$) **group** collapses it to the
single flavour-blind number $N_f$ — which can carry only *one* dressing. So to give components
different dressings you either avoid the group (split into separate diagrams) or keep the group
but make its $\delta$ dress selected components with distinctly-named dressings.

## Per-flavour split

Here u and d are two **separate fields**, each its own diagram with its own dressing. There is
no SU($N_f$) flavour group at all — only colour SU(3) is folded. The generator
`tests/gen/gen_flavour_split.wls` builds a scalar (sigma) two-point self-energy from a quark
loop, summed over the two flavours, with hand-rolled Yukawa rules:

```mathematica
(* one scalar sigma + TWO quark flavours, each a field of its own *)
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qub[p, {di, ci}], qu[p, {di, ci}]},
                  {qdb[p, {di, ci}], qd[p, {di, ci}]}}|>;

(* the only difference between u and d is the dressing name on the propagator / Rdot insertion *)
dressed = (expr /. myRules) //. {
  dressing[GammaN, {qub, qu}, 1, {m1_, m2_}] :> Gu[Sqrt[sp[m2, m2]]],
  dressing[GammaN, {qdb, qd}, 1, {m1_, m2_}] :> Gd[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qub, qu}, 1, {m1_, m2_}] :> GuDot[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qdb, qd}, 1, {m1_, m2_}] :> GdDot[Sqrt[sp[m2, m2]]]};

SetNc[3];                          (* colour SU(3); no flavour group — no flavour heads in the net *)
net = FromFunKit[dressed];
ntk = NumTrace[net, "Frame" -> frame, "Args" -> {l1, cos1, p}];
MakeNTKernel[ntk, "gen_flavour_split_num.cpp",
  "Flavour_split_num_kernel.hh", "Flavour_split_num_kernels.hh",
  "Dressings" -> {Gu, Gd, GuDot, GdDot}];        (* four INDEPENDENT runtime dressings *)
```

The kernel returns `(u-loop with Gu) + (d-loop with Gd)` and exposes `Gu, Gd, GuDot, GdDot`
as independent runtime parameters. `compare_flavour_split.cpp` (ctest `flow_flavour_split`)
validates it against the FORM oracle.

## Group-diagonal dressing

The same physics, but built as a **single** quark loop carrying an SU(2) flavour index whose
propagators carry a *flavour-diagonal* dressing $\mathrm{diag}(G_u, G_d)$ instead of a
flavour-blind $\delta$. The new heads are

```mathematica
ntSUNDiagFund[N, i, j, spec, s]   (* fundamental δ_ij dressing selected components *)
ntSUNDiagAdj [N, a, b, spec, s]   (* adjoint     δ^ab dressing selected components *)
```

`spec` is a **rules list** `{c1 -> name1, c2 -> name2, …, Default -> defName}`: the `ci` are
**1-based** component indices (1..$N$ fundamental, 1..$N^2-1$ adjoint) and the `namei` are
distinctly-named **scalar** dressing symbols evaluated at the scale `s`. A `Default -> defName`
rule dresses every unnamed component; components with neither a rule nor a `Default` are
**dropped** (contribute nothing — no dead terms). The SU($N$) trace then folds — via
`sun_value_dressed` (`network/sun_net.hpp`) — to $\sum_a c_a\, Z_a(s)$ over those named runtime
scalar dressings (each an ordinary kernel parameter), **and the Dirac/colour trace is computed
once**, not once per component. From `tests/gen/gen_flavour_ingroup.wls`:

```mathematica
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qb[p, {di, ci, fi}], q[p, {di, ci, fi}]}}|>;   (* ONE flavour doublet *)

(* u = component 1, d = component 2 — both named, nothing drops *)
diagFlav[fb_, f_, nmU_, nmD_, m_] := NumTracer`ntSUNDiagFund[2, fb, f, {1 -> nmU, 2 -> nmD}, Sqrt[sp[m, m]]];
myRules = {
  (* Yukawa vertex: flavour-BLIND δ (deltaFundFlav -> ntSUNDeltaFund) *)
  GammaN[{sig, qb, q}, …] :> I gamma5[dqb, dq] deltaFundCol[cqb, cq] ntSUNDeltaFund[2, fqb, fq],
  (* propagator: flavour-DIAGONAL dressing diag(Gu, Gd) in place of the blind δ × dressing *)
  Propagator[{q, qb} | {qb, q}, …] :> deltaFundCol[c1, c2] diagFlav[f1, f2, Gu, Gd, m2] (-I) gamma[mu, d1, d2] vec[m2, mu],
  Rdot[{qb, q}, …]               :> deltaFundCol[c1, c2] diagFlav[f1, f2, GuDot, GdDot, m2] (-I) gamma[mu, d1, d2] vec[m2, mu]};

MakeNTKernel[ntk, …, "Dressings" -> {Gu, Gd, GuDot, GdDot}];   (* ordinary scalar dressings *)
```

This is the in-group counterpart of the split flow and **reproduces it exactly** — both kernels
now share the same `(Gu, Gd, GuDot, GdDot)` scalar signature:
`compare_flavour_ingroup.cpp` (ctest `flow_flavour_ingroup`) checks the two kernels agree
pointwise.

### The adjoint case: a gluon condensate

`ntSUNDiagAdj` does the same for an adjoint index — the natural home of a condensed gluon
dressing $Z_A^a$. A condensate lives in the **Cartan** directions ($\lambda_3, \lambda_8$ for
SU(3)), so only components 3 and 8 need distinct dressings and the other six drop out:

```mathematica
ntSUNDiagAdj[3, a, b, {3 -> ZA3, 8 -> ZA8}, scale]   (* dress λ3, λ8; drop the other 6 *)
```

`tests/gen/gen_gluon_condensate.wls` emits three kernels: `_cartan` (the above — only 3 & 8
survive), `_diag` (`{3 -> ZA3, 8 -> ZA8, Default -> ZAoff}`, the other 6 collapsed into a shared
`ZAoff`), and the colour-blind `_plain` (`ntSUNDeltaAdj[3, a, b]` × a scalar). The test
(`flow_gluon_condensate`) validates two properties: the **collapse** — feeding `_diag` one
function for `ZA3 = ZA8 = ZAoff` reproduces `_plain` exactly ($\sum_a Z^2 = 8\,Z^2$) — and the
**drop** — `_cartan` equals `_diag` with `ZAoff ≡ 0`, i.e. dropping a component and dressing it
with a zero function agree.

```{note}
Group-diagonal dressings leave `sun_value_cx` (the plain colour fold) untouched, so any flow
that does *not* use them regenerates byte-identical.
```

## Verify

```bash
cd numtracer
# regenerate (needs wolframscript + FunKit + FORM):
wolfram -script tests/gen/gen_flavour_split.wls
wolfram -script tests/gen/gen_flavour_ingroup.wls
wolfram -script tests/gen/gen_gluon_condensate.wls
# validate the committed kernels (no Wolfram needed):
cmake --build build -j4
ctest -R "flow_flavour_split|flow_flavour_ingroup|flow_gluon_condensate" --output-on-failure
```

For the design — how the colour fold keeps a group-diagonal $\delta$ as a polynomial and how a
diagram that disconnects into several traces is emitted as a product — see
[the numeric engine](../internals/numeric-engine.md) and [front-end & codegen](../internals/codegen.md).
