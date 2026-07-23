# Dressed propagator numerators

Tutorial 6 dressed a colour/flavour **index** — a group-diagonal $\delta$ that folds to a
[`SUNPoly`](dressed-flavour.md). This tutorial dresses the **Dirac structure** of a propagator
numerator, which folds to the analogous [`DPoly`](../internals/numeric-engine.md). It is the last
of the three "atom-carrying polynomial" layers ([`MPoly`](../internals/numeric-engine.md) for
$1/k^2$, `SUNPoly` for colour, `DPoly` for Dirac numerators) and the one closest to everyday fRG:
a dressed quark propagator.

## The problem: dressed numerators multiply out

A dressed quark propagator numerator is not a single Dirac structure but a **sum** of them, each
with its own *runtime* dressing:

$$
S(p) \;=\; M_q\,\mathbb{1} \;+\; Z(p)\,\slashed p ,
$$

a mass term plus a wave-function-dressed slash. A diagram with $D$ such numerators, if you
distribute the sums, is $2^D$ separate traces — one per choice of structure in each numerator.
That is the blow-up NumTracer avoids: it keeps each numerator **eager** (a *slot* holding the
numerator's list of structure options), contracts the chain **once**, and collects the result
into a `DPoly` — a polynomial whose variables are the opaque dressing calls and whose coefficients
are the kinematic `MPoly`s the engine already computes. The dressings never enter the trace
arithmetic; they ride along as atom-ids, exactly like `MPoly`'s $1/k^2$ atoms. So the Dirac and
Lorentz work happens a single time regardless of how many structures each numerator carries.

## The runnable program

We take a quark bubble with **two** dressed numerators, closed by the gluon metric:

$$
\operatorname{tr}\!\big(\gamma^\mu\, S(p)\, \gamma^\nu\, S(q)\big)\,\delta_{\mu\nu},
\qquad S(p)=M_q\mathbb{1}+Z(p)\slashed p,\quad S(q)=M_q\mathbb{1}+Z(q)\slashed q .
$$

Distributed, that is $2\times2=4$ traces; collected, it is **one** `DPoly`. The program builds
the collected `DPoly` once and checks it against the explicit distributed sum at random kinematics
and random dressing values.

```cpp
#include <numtracer.hpp> // numeric_value_dressed / DSlot / DPoly / eval + network::dgamma / dslash

using numtracer::Cx;
namespace nm = numtracer::numeric;  // MPoly / DPoly / DSlot / numeric_value / numeric_value_dressed
namespace net = numtracer::network; // DiracNet / dgamma / dslash / NetVal builders

// Frame symbols: momenta p (symbols 0..3) and q (symbols 4..7); comp[vid] = its 4 MPoly components.
const int nsym = 8;
nm::LorentzEnv env(nsym); // bind the 8-symbol space once; MPoly/DPoly are minted from it
std::vector<std::array<nm::MPoly, 4>> comp(2);
for (int mu = 0; mu < 4; ++mu) {
  comp[0][mu] = env.var(mu);     // p_mu
  comp[1][mu] = env.var(4 + mu); // q_mu
}

// Lorentz half: the two free gluon legs (ids 100, 101) meet through one metric.
nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};

// The dressed numerators. A DSlot is a numerator's list of structure options; a DSlotOpt is
//   { coeff , dressing-atom ids , Dirac-token chain , Lorentz-net factors }.
//   S(p): option 0 = 𝟙 dressed by atom 0 ("Mq");   option 1 = p̸ dressed by atom 1 ("Z(p)").
//   S(q): option 0 = 𝟙 dressed by atom 0 ("Mq");   option 1 = q̸ dressed by atom 2 ("Z(q)").
// An empty token chain IS the spinor identity 𝟙; the net-factor list is empty here (see below).
// Atom 0 is SHARED: both mass terms are the same runtime Mq.
nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}},
                nm::DSlotOpt{Cx{1, 0}, {1}, {net::dslash({{1.0, 0}})}, {}}}; // p̸ = 1·comp[0]
nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}},
                nm::DSlotOpt{Cx{1, 0}, {2}, {net::dslash({{1.0, 1}})}, {}}}; // q̸ = 1·comp[1]

// The dressed Dirac chain in trace order: dtfix = a fixed factor, dtslot(i) = the i-th slot.
std::vector<nm::DChainTok> dchain = {nm::dtfix(net::dgamma(100)), nm::dtslot(0),
                                     nm::dtfix(net::dgamma(101)), nm::dtslot(1)};

// Collect: ONE contraction, no 2^D blowup.
nm::DPoly dp = env.numeric_value_dressed(dchain, {sP, sQ}, lor, comp, /*atomDen*/ {});
```

Evaluating the `DPoly` (`eval(dp, x, atomVal, drVal)` — kinematics `x`, $1/k^2$ atom values
`atomVal`, dressing values `drVal`) against the distributed reference over 5000 random points
prints:

```text
Dressed quark bubble  tr( γ^μ S(p) γ^ν S(q) ) δ_{μν},  S = Mq·𝟙 + Z·slash
  DPoly dressing monomials (one trace collected) : 2   (distributed would be 4 traces)
  collected == distributed over 5000 random points : worst |Δ| = 0.00e+00  (0 bad)
ALL TESTS PASSED
```

## Reading the result

**Why two monomials, not four?** The four structure choices carry dressing products
$M_q^2,\; M_q Z(q),\; Z(p)M_q,\; Z(p)Z(q)$. But the two mixed choices are traces of an **odd**
number of $\gamma$'s — $\operatorname{tr}(\gamma^\mu\gamma^\nu\slashed q)$ and
$\operatorname{tr}(\gamma^\mu\slashed p\,\gamma^\nu)$ — which vanish identically. So collection
keeps only $M_q^2$ (from $\mathbb{1},\mathbb{1}$) and $Z(p)Z(q)$ (from $\slashed p,\slashed q$),
and the two structurally-zero traces are dropped **for free**. The `DPoly` is
$M_q^2\,\big[\operatorname{tr}(\gamma^\mu\gamma^\nu)\delta_{\mu\nu}\big] + Z(p)Z(q)\,\big[\operatorname{tr}(\gamma^\mu\slashed p\,\gamma^\nu\slashed q)\delta_{\mu\nu}\big]$,
each bracket a kinematic `MPoly` computed once.

**The data model.** A `DPoly` is a sorted list of `(dressing monomial, MPoly)` pairs. A **dressing
monomial** (`DMono`) is a sorted multiset of dressing-atom ids — here `{0,0}` for $M_q^2$ and
`{1,2}` for $Z(p)Z(q)$ — merged on multiplication exactly like `MPoly`'s inverse-atom multiset.
The coefficient of each dressing monomial **is a full `MPoly`**, so `DPoly` reuses
`MPoly::operator*`/`operator+` verbatim and only diagrams that actually carry a dressed structure
sum pay for the dressing layer. See
[the numeric engine](../internals/numeric-engine.md) for how the id-namespaces relate.

**Beyond propagator numerators.** A `DSlotOpt`'s structure is described by two lists, because a
collected slot is not restricted to the $\{\mathbb{1}, \slashed p\}$ pair above:

* `toks` — a **Dirac-token chain** spliced in place of the slot (spinor `din → dout`). Any free
  Lorentz index it carries (a `dgamma(μ)`, an open `dcomm` leg) is an *open leg*.
* `netFacs` — **Lorentz-net factors** (a vector $p^\mu$, a metric $g^{\mu\nu}$, …) carrying any
  remaining open legs, appended to the surrounding net.

Every option in one slot shares the same open-leg set, so the surrounding net closes a fixed set of
legs whatever structure is chosen. That one form covers any leg count $k \ge 0$: $k=0$ is the
propagator numerator of this tutorial (`toks = {}` for $\mathbb{1}$, `{dslash(p)}` for
$\slashed p$, both with `netFacs = {}`); $k=1$ is a quark–gluon vertex ($\gamma^\mu$,
$\slashed p_1\gamma^\mu$, $\sigma^{\mu\nu}\slashed p_\nu$, or $p^\mu\cdot\mathbb{1}$ with the leg on
a `netFacs` vector); $k=2$ a two-gluon vertex, and so on. `tests/test_dpoly.cpp` cases G–J walk that
ladder against explicit distributed sums.

```{admonition} How the front-end produces slots
:class: note
You rarely hand-build `DSlot`s: the Mathematica front-end keeps a dressed numerator eager and
emits the slot chain, and the generated kernel lowers each dressing atom to a `dress` env leaf
(an opaque runtime value the kernel evaluates once, like an `inv` $1/k^2$ leaf). The committed
`Zq` collection flow (`tests/gen/gen_zq_collect.wls`, validated by `compare_zq_collect.cpp`) is
the full-pipeline version of this tutorial; `tests/test_dpoly.cpp` exercises the eager collection
directly on the harder $\sigma$-vertex chains.
```

## Verify

```bash
cd Tutorials
cmake -S . -B build            # configure once
cmake --build build --target dressed_numerators
./build/dressed_numerators     # self-checks collected == distributed; exits non-zero on mismatch
```

The engine-level correctness gates for the whole dressing-collection layer live in the
`numtracer/` build: `ctest -R dpoly` (arithmetic, the no-dressing regression guarantee, and
dressed contractions up to the 147-vertex $\sigma$ cases) and `ctest -R flow_zq_collect` (the
committed full-pipeline `Zq` collection kernel against the FORM oracle).
