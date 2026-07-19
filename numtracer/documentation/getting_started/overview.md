# Overview

NumTracer assembles and evaluates the tensor networks that appear in quantum-field-theory
loop integrands, and **generates the scalar kernels** that a numerical integrator then runs.
It does the tensor algebra in plain compile-time C++ — no symbolic-algebra runtime, no external
tools — and emits a flat, straight-line real-arithmetic kernel that a consumer compiles directly.

## What it computes

A loop integrand is a product of tensors — gamma matrices, momentum vectors, metrics,
projectors, colour generators — contracted together (summed over shared indices) down to a
single number that still depends on a few runtime quantities: momentum magnitudes, angles, and
propagator dressings. That number is fed to an integrator and evaluated at hundreds of
thousands of grid points.

NumTracer is a *general* tensor-tracing engine. The physics is entirely in the network you
hand it; the engine only contracts indices and folds the result. The numeric path handles
Lorentz, Dirac, and SU($N$) colour *and* flavour structure, dresses individual flavours or
colour/flavour components differently (per-flavour split and group-diagonal dressings), and factors
a diagram that disconnects into several closed traces into a product — see
[Key concepts](concepts.md) and the [dressed-flavour tutorial](../tutorials/dressed-flavour.md).

```{admonition} Is this for me?
:class: tip
If you have a tensor network built from metrics, vectors, projectors, gamma matrices, and/or
SU($N$) factors — in a **4D Euclidean** setting — NumTracer will contract it to a scalar kernel,
whatever theory it comes from. The worked flows here are fRG flows for Yang–Mills and QCD, but that
is just the authors' application. Before assuming it is (or isn't) for you, read
[Scope & conventions](scope-and-conventions.md) for the exact boundary, then
[Bring your own network](bring-your-own-network.md) for a domain-neutral hello-world and the
dictionary that maps your objects onto the engine's.
```

## The numeric path, end to end

NumTracer generates kernels by **numeric tracing** — a build-time generator driven from a small
Mathematica front-end:

1. **Front-end.** `NumTrace` analyses a tensor network written in DSL heads (`ntMetric`,
   `ntVec`, `ntTransProj`, `ntSUNf`, `ntSUNDeltaAdj`, `ntDress`, …) into a list of diagrams — each a
   scalar coefficient times a contraction — together with the loop frame that fixes every
   momentum's components and the env-id layout. `FromFunKit` imports a FunKit flow into the
   same DSL. Vertex structure-sums are kept eager, so no product-of-sums monomial blow-up ever
   forms; the same holds for dressed propagator numerators, which are collected into a single
   trace rather than distributed into `2^D` diagrams (`"DressingCollection"`, on by default).

2. **Numeric contraction.** `MakeNTKernel` emits a small C++ generator and runs it. The
   generator contracts each diagram numerically over the frame: the Dirac trace as 4×4 chiral
   matrix products (including γ5; no `(2n−1)!!` Wick blow-up), the Lorentz network by bounded
   index elimination (no `2^np` projector-mask blow-up), and the colour factor folded to a
   number. Each diagram becomes a small polynomial in the frame's scalar symbols.

3. **Lowering.** Common-subexpression elimination and Horner factoring lower that polynomial
   to a flat straight-line kernel — trace functions over the frame symbols, a `fill` that
   computes the symbols once per call, and the per-diagram assembly. The committed kernel is
   plain C++ that depends only on a tiny runtime-support header.

The same contraction primitives are exposed directly in C++, both for hand-built traces and
as the numeric oracle the generated kernels are validated against.

## Why generate, rather than evaluate symbolically

Doing the full tensor contraction inside the consumer's compiler is correct but expensive: a
trace of several transverse projectors expands explosively in the frame-component basis, and the
compiler never reclaims intermediate memory. Running the *same* contraction as a build-time step —
numerically, over a fixed frame — sidesteps both: the generator runs in seconds and tens of
megabytes, and the consumer only ever compiles the small, flat result.

```{admonition} Relationship to symbolic tracers
:class: note
This is the one place we compare NumTracer to a symbolic tensor-algebra system such as FORM (used,
via FormTracer, only as a validation oracle for the test suite). Such a tool does the tensor
algebra symbolically ahead of time and emits a flat polynomial in the scalar products; NumTracer
produces the *same kind* of flat scalar kernel, but does the contraction numerically in C++ over a
fixed frame. In practice it generates a kernel **~80–175× faster**, and the generated kernel is
**competitive with or faster than** the symbolic one: on the quark–gluon vertex `ZAqbq{1,4,7}_147`
it runs at **0.96× / 0.99× / 0.62×** the FORM kernel's time (`tests/refshim/bench_aqbq147.cpp`),
and on the pure-gauge `ZA3_147` at 1.01×. See [PERFORMANCE.md](../../PERFORMANCE.md) for the
per-flow table. Nothing downstream depends on that tool: the emitted kernel is self-contained C++.

A fixed-frame contraction was long assumed to be structurally unable to do the partial-fractioning
(integration-by-parts-like) step a symbolic tracer performs on scalar products before the frame is
substituted, and that was believed to leave an irreducible residual. It does not: the cancellation
can be done *in* the frame, by exact polynomial division of a shifted-line propagator denominator
into the numerator that contains it (`divThroughPolyAtoms`, see
[the numeric engine](../internals/numeric-engine.md)). That closed the residual, and improved
accuracy at the same time — the division and the terms that cancel against it both disappear.
```

Next, [Key concepts](concepts.md) builds the mental model the
[tutorials](../tutorials/index.md) rely on. See the [C++ API](../doxygen/NumTracer/html/index)
for the full reference.
