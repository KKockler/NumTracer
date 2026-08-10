# step-10: From the Wetterich equation to a network

*Builds on: [step-07](step-07.md) · Built on by: [step-12](step-12.md) · Tags: `funkit`, `frg` ·
**Tier C** (the full fRG toolchain)*

````{admonition} What this step needs
:class: important
[FunKit](https://github.com/satfra/FunKit), TensorBases, DiFfRG **and** FORM, in addition to a
Wolfram kernel. This is the first step that needs the whole toolchain, and it is therefore **not**
wired into `ctest` — the script is here to be read, and to be run if you have the pieces:

```bash
NT_TUTORIAL_OUT=/tmp/out wolframscript -script Tutorials/step-10-wetterich/quark_selfenergy.wl
```

It writes only to `NT_TUTORIAL_OUT` and modifies nothing in the repository.

FORM appears here because TensorBases uses it to compute the tensor *bases* — never to do
NumTracer's trace.
````

## Introduction

Every step so far handed `NumTrace` a network we wrote out by hand. In a real calculation nobody
does that. The network is **derived**, and this step shows the derivation for the running
example — the quark two-point flow, whose numerator [step-05](step-05.md) assembled by hand.

### Where a flow comes from

The Wetterich equation is an exact identity for the scale derivative of the effective average
action $\Gamma_k$:

$$
\partial_t \Gamma_k \;=\; \tfrac12\,\mathrm{STr}\Big[\big(\Gamma_k^{(2)} + R_k\big)^{-1}\,\partial_t R_k\Big].
$$

The supertrace runs over all fields; $R_k$ is the regulator that suppresses modes below the scale
$k$. To get the flow of a particular vertex you differentiate both sides with respect to the
corresponding fields. Two derivatives $\delta^2/\delta\bar q\,\delta q$ give the flow of the quark
two-point function.

That differentiation is mechanical and is what FunKit does. The result is a sum of diagrams, each a
product of propagators, vertices and one regulator insertion, with loop momentum routed and all
indices carried symbolically. NumTracer's job begins only after that.

### The five stages

```text
WetterichEquation
      │  FTakeDerivatives · FTruncate · FRoute        (1) differentiate
      ▼
  a sum of diagrams, indices open
      │  FTerm[projector] ** …                        (2) project to a scalar
      ▼
  a scalar flow, abstract dressing heads
      │  dressingRules · PropParam                    (3) name dressings, reduce scalars
      ▼
  a traced flow in FunKit's language
      │  FromFunKit                                   (4) import
      ▼
  a NumTracer DSL network            ← steps 06-09 started HERE
      │  NumTrace · MakeNTKernel                      (5) trace and emit
      ▼
  a C++ kernel
```

Stages 1–4 are FunKit and TensorBases; only stage 5 is NumTracer. It is worth being clear about
that boundary, because most "NumTracer gives the wrong answer" reports turn out to live in stages
2 or 3.

## The commented program

`Tutorials/step-10-wetterich/quark_selfenergy.wl`.

### Loading the toolchain

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
:start-after: "@snip begin: load"
:end-before: "@snip end: load"
```

`qcd_setup.m` is **not part of NumTracer**. It is the per-theory model file: field content,
truncation, which tensor basis each vertex uses, and the dressing rules. In your own project this
is your `QCD.wl`. The copy under `numtracer/tests/` is a trimmed version of the production
`QCD_Nf2` setup and is worth reading alongside this script.

```{admonition} `SetDirectory[$repo]` is load-bearing
:class: note
TensorBases resolves its basis cache (`TBCache/`) **relative to the working directory**. Run the
script from elsewhere and it will recompute bases from scratch through FORM — minutes instead of
seconds — or fail to find them at all.
```

### (1) Differentiate

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
:start-after: "@snip begin: derive"
:end-before: "@snip end: derive"
```

Three operations, and the order matters:

* **`FTakeDerivatives`** does the functional differentiation. The external index labels `i1`, `i2`
  are placeholders; the association `fRGqbq["1-Loop"]["ExternalIndices"]` records what they were
  bound to, which the projector needs in stage 2.
* **`FTruncate`** discards everything outside the truncation declared in the setup. A vertex you
  did not declare simply does not appear — silently, by design.
* **`FRoute`** assigns loop-momentum routing: which internal line carries $l$ and which carries
  $l-p$. This is **leg-order sensitive**, which is why `qcd_setup.m` pins FunKit's backend before
  any derivative is taken. Change the routing and every frozen test integrand changes — the same
  integral, but a differently-shaped integrand.

### (2) Project

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
:start-after: "@snip begin: project"
:end-before: "@snip end: project"
```

The flow of a two-point function is a *matrix* in Dirac and colour space. To extract a scalar you
contract with the dual of the basis element you want — element 1 of the `"qbq"` basis, the quark
wave function.

Three details, each of which has cost somebody a week:

* **`**` is FunKit's non-commutative product.** Ordinary `*` would commute the projector into the
  middle of a Dirac chain.
* **The index order is `{i2, i1}` — reversed** relative to the derivative order `{qb[i1], q[i2]}`.
* **The leading minus** compensates the sign convention of the `qbq` dressing: `qcd_setup.m` maps
  `dressing[GammaN, {qb, q}, 1, …] :> -Zq[…]`. Drop it and the whole flow changes sign.

None of these three is checkable by the engine. [step-11](step-11.md) is about how to *assert* them
instead of hoping.

### (3) Dress and reduce

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
:start-after: "@snip begin: dress"
:end-before: "@snip end: dress"
```

Two substitution layers, in this order.

**`dressingRules`** maps the abstract heads the derivative produced —
`dressing[GammaN, {A, A}, 1, {p1, p2}]` — onto the named functions the kernel will take as
parameters, and inserts the regulator. From `qcd_setup.m`:

```mathematica
dressing[GammaN, {A, A},  1, {p1_, p2_}] :>  ZA[Sqrt[sp[p2, p2]]] sp[p2, p2],
dressing[GammaN, {qb, q}, 1, {p1_, p2_}] :> -Zq[Sqrt[sp[p2, p2]]],
dressing[InverseProp, {A, A}, 1, {p1_, p2_}] :>
    ZA[Sqrt[sp[p2, p2]]] sp[p2, p2] + RB[k^2, sp[p2, p2]] ZA[evP],
```

````{admonition} An unmatched dressing rule does not fail — it leaks
:class: warning
A `dressing[…]` head that no rule matches is **not** an error. It survives the substitution, rides
through `FromFunKit` as an opaque scalar, and is emitted into the kernel as a call to an undefined
function:

```cpp
const auto _interp7 = dressing(GammaN, List(qb, q), 1, ...);   // does not compile
```

The production `QCD_Nf2` code documents the classic instance: FunKit's `FOrderFields` canonicalises
a four-fermi vertex to `{qb, qb, q, q}`, so a rule written for the `{qb, q, qb, q}` spelling **never
fires**. A `dressing(...)` in a generated file always means a rule that did not match — see
[step-22](step-22.md).
````

**`PropParam`** does the scalar kinematic reduction: `sp[l1, p1] -> l1 p cos1` and friends. This is
the **scalar half only**. The tensor half gets its kinematics from the frame, and the two must
agree — which for multi-angle flows is what `"AngleDefs"` is for ([step-07](step-07.md),
[step-14](step-14.md)).

### (4) Import, and (5) trace

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
:start-after: "@snip begin: import"
:end-before: "@snip end: import"
```

`FromFunKit` rewrites the traced flow into the DSL of [step-08](step-08.md): gammas become
`ntGamma`, momentum legs `ntVec`, scalar products `ntSP`, and SU($N$) tokens get their rank baked
in from `Nc` and `Nf`. **From here on it is exactly the object steps 06–09 built by hand.**

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
:start-after: "@snip begin: emit"
:end-before: "@snip end: emit"
```

The diagram-count assertion is worth keeping. The quark two-point flow has exactly **two** diagrams
— the regulator insertion sits on the gluon line in one and on the quark line in the other. Getting
1, or 5, means the truncation or the routing is not what you think.

## Results

```text
step-10: 1-loop external indices = {i1 -> {p1, {d1, A1, F1}}, i2 -> {-p1, {d2, A2, F2}}}
step-10: diagrams = 2   (expected 2: the regulator insertion sits on the gluon line in one
                         and on the quark line in the other)
wrote header: .../Zq_tutorial_kernels.hh (2716 bytes)
wrote kernel: .../Zq_tutorial_kernel.hh
```

The external-index line is the first thing to read in any flow. It says $i_1$ carries momentum
$p_1$ with indices $(d_1, A_1, F_1)$ — Dirac, colour, flavour — and $i_2$ carries $-p_1$. That
$-p_1$ is momentum conservation, and it is the reason the projector in stage 2 needed the reversed
index order.

### The generated kernel

```cpp
static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k,
                          const auto& ZA, const auto& Zq, const auto& Mq,
                          const auto& ZAqbq1, const auto& dtZA, const auto& dtZq)
{
  double fenv[nenv];
  const double dr_0 = Mq(l1);
  const double dr_1 = -sqrt(powr<-1>(powr<2>(l1))) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
  const double dr_2 = Mq(sqrt(powr<2>(l1) - 2.*cos1*l1*p + powr<2>(p)));
  const double dr_3 = …;
  fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3);
  const auto _interp1 = ntRe(tr0(fenv));
  …
}
```

Compare this with the toy kernel of [step-06](step-06.md); the shape is identical and the content is
now real physics.

* **The signature is `"Args"` then `"Dressings"`**, in order. Six dressings, each an opaque callable.
* **`dr_0 … dr_3` are *dressed* quantities entering the trace itself.** This is new. A dressed quark
  propagator numerator is $M_q\mathbb{1} + Z(p)\slashed p$ — a sum of Dirac structures with
  *runtime* coefficients — so those coefficients have to reach the trace. They are passed into
  `fill` as extra env slots rather than distributing the diagram into $2^D$ copies. That mechanism
  is [step-18](step-18.md); here you are seeing its footprint.
* **`sqrt(powr<2>(l1) - 2.*cos1*l1*p + powr<2>(p))` is $|l-p|$**, written out in the frame's
  scalars. The scalar reduction of stage 3 produced it.
* **`nenv = 7` and one trace function.** The whole two-diagram flow collapsed to a single
  `tr0(fenv)`. The dressings did not multiply the number of traces.
* **`ntRe(...)`** — the trace is complex-valued in general (the $i$'s in the gammas and the colour
  factors), and the physical integrand is its real part. [step-21](step-21.md) explains how the
  generator decides that this is legitimate.

## Possibilities for extensions

1. **Read the routing.** Print `fRGqbq["1-Loop"]["Expression"]` before the projection and identify
   the two diagrams. Which propagator carries $l$ and which $l-p$? Then flip the backend pin in
   `backend_pin.m` and look again — same integral, different integrand.

2. **Drop the minus sign.** Remove the `-` in the `FTerm[…]` and regenerate. The kernel is
   perfectly valid and the flow runs the wrong way. Nothing warns you. This is the single most
   common sign error in the whole workflow.

3. **Force a leaked dressing.** Comment out the `dressing[GammaN, {A, qb, q}, 1, …]` rule in a
   *copy* of `qcd_setup.m` and regenerate. Find the `dressing(...)` call in the emitted file and
   note how far downstream the failure appears relative to where the mistake was.

4. **Truncate harder.** Add `//. {ZAqbq1[__] :> 0}` after `dressingRules` and watch the diagram
   count and kernel size drop. This is exactly the `dropAqbq47` / `dropFourQuark` idiom the
   production `QCD_Nf2` code uses to restrict a trace ([step-14](step-14.md)).

5. **Compare against the committed kernel.** `numtracer/tests/gen/Zq_num_kernel.hh` is the same
   flow, generated against the test shim. Diff the two traces headers. They should agree up to the
   namespace and the support API — and if they do not, one of you has changed the physics.

## The plain program

```{literalinclude} ../../../Tutorials/step-10-wetterich/quark_selfenergy.wl
:language: mathematica
```

Next, [step-11](step-11.md): the projector in stage 2 — where it comes from, what its index order
means, and how to assert that it is the one you meant.
