# Per-flavour and per-component dressings

This is the most physics-heavy tutorial; the [Glossary](../getting_started/glossary.md) defines
*dressing*, *flavour*, *isospin*, and *condensate*. Its full-flow artifacts are the Mathematica
generators (`tests/gen/gen_flavour_split.wls`, `gen_flavour_ingroup.wls`, `gen_gluon_condensate.wls`)
and their committed kernels — validated by the ctests named below. The **group-diagonal fold** at
the heart of the in-group approach is also directly callable from standalone C++, with no
Mathematica involved — the [pure-C++ section](#in-pure-c-the-group-diagonal-fold) below is a
runnable worked example.

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

## In pure C++: the group-diagonal fold

The Mathematica above is only front-end sugar. The mechanism it drives — the SU($N$) trace folding
to $\sum_a c_a\,D_a$ over named dressings instead of to one flavour-blind number — is the numeric
engine's colour/flavour fold (`network/sun_net.hpp`), callable directly. There are two entry points:

- `sun_value_cx(net)` — the **plain** fold: a fully-contracted SU($N$) net → **one number** `Cx`.
- `sun_value_dressed(net)` — the **group-diagonal** fold: if the net carries a `diagFund`/`diagAdj`
  factor it returns a **`SUNPoly`** — a small polynomial $\sum_t c_t \prod D^{\mathrm{dr}}$, each
  term a `SUNTerm{coeff, dress}` (`dress` = the list of dressing-ids in that monomial). A net with
  no diagonal factor comes back as a single constant term, exactly `sun_value_cx`.

A `diag` factor is a $\delta$ tagged with a `comp2dr` map — `comp2dr[v]` is the dressing-id for
component `v`, or `-1` to **drop** it. This is the C++ image of the Mathematica `spec` rules list.
The runnable program `Tutorials/06-dressed-flavour/dressed_flavour.cpp` builds the flavour story of
this tutorial (u/d in an SU(2) doublet) at the engine level and self-checks it — its core:

```cpp
#include <numtracer.hpp> // sun_value_dressed + SUN::diagFund / diagAdj / deltaFund / deltaAdj

using namespace numtracer;          // Cx, approx
namespace net = numtracer::network; // SUNPoly / SUNTerm / sun_value_dressed / SUN

enum { i, j }; // the two ends of the closed δ loop

// Evaluate a SUNPoly at a dressing-id -> value assignment D: Σ_t coeff_t · Π_{id∈dress} D(id).
static Cx eval_poly(const net::SUNPoly &p, double (*D)(int)) {
  Cx s{0, 0};
  for (const net::SUNTerm &t : p) {
    Cx c = t.coeff;
    for (int id : t.dress) c = c * Cx{D(id), 0.0};
    s = s + c;
  }
  return s;
}

int main() {
  // A flavour δ-loop over an SU(2) doublet, its δ made flavour-DIAGONAL: component 0 (u) carries
  // dressing-id 0, component 1 (d) carries id 1. (comp2dr is 0-based: {0, 1} = "u->id0, d->id1".)
  // sun_value_dressed folds the closed loop diagFund·deltaFund to the SUNPoly  D_u + D_d.
  const net::SUNPoly ud = net::sun_value_dressed(
      {net::SUN::diagFund(2, i, j, {0, 1}), net::SUN::deltaFund(2, j, i)});

  auto broken = [](int id) { return id == 0 ? 2.0 : 5.0; }; // D_u = 2, D_d = 5 (broken isospin)
  auto ones = [](int) { return 1.0; };                      // all dressings equal to 1
  std::printf("   D_u + D_d           = %g   (D_u=2, D_d=5 -> 7)\n", eval_poly(ud, broken).re);
  std::printf("   collapse D_u=D_d=1  = %g   (-> flavour-blind N_f = 2)\n", eval_poly(ud, ones).re);

  // DROP: the same loop dressing ONLY component 0 (id -1 drops component 1) folds to just D_u.
  const net::SUNPoly u_only = net::sun_value_dressed(
      {net::SUN::diagFund(2, i, j, {0, -1}), net::SUN::deltaFund(2, j, i)});
  std::printf("   drop d (only u)     = %g   (-> D_u = 2)\n", eval_poly(u_only, broken).re);
  return 0; // part B (the adjoint gluon condensate) follows in the source
}
```

Building the tutorial (`cmake -S Tutorials -B Tutorials/build && cmake --build Tutorials/build`) and
running `./Tutorials/build/dressed_flavour` prints:

```text
A. fundamental u/d doublet (SU(2) flavour)
   D_u + D_d           = 7   (D_u=2, D_d=5 -> 7)
   collapse D_u=D_d=1  = 2   (-> flavour-blind N_f = 2)
   drop d (only u)     = 2   (-> D_u = 2)
```

Those three lines *are* the physics: `D_u + D_d` is the broken-isospin loop (two flavours dressed
independently), the collapse `D_u = D_d = 1 → N_f` recovers the flavour-blind number a plain
`sun_value_cx` would give (this is what `flow_gluon_condensate`'s collapse property checks in the
adjoint), and the drop reproduces a condensate living on only some components. The Dirac/Lorentz
trace multiplying this fold is computed **once** — the `SUNPoly` supplies the per-flavour
coefficients as a polynomial, not one diagram per component. The tutorial's **part B** does the
adjoint case identically with `SUN::diagAdj(N, a, b, comp2dr)` (dimension $N^2-1$): dressing only
components 3 and 8 with `-1` elsewhere is the Cartan gluon condensate `{3,8}` of the Mathematica
example (`full loop, all Z=1 = 8`, `Cartan {3,8} = 2`).

## Verify

```bash
cd numtracer
# the pure-C++ fold above (no Wolfram needed):
cmake --build build -j4 && ctest -R col_numeric --output-on-failure
# regenerate the full flows (needs wolframscript + FunKit + FORM):
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
