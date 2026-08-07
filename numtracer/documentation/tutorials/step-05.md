# step-05: A full diagram, and what lowering does to it

*Builds on: [step-02](step-02.md), [step-03](step-03.md), [step-04](step-04.md) · Built on by:
[step-06](step-06.md) · Tags: `capstone`, `codegen`, `cse` · **Tier A** (a C++20 compiler, nothing
else)*

## Introduction

This is the capstone of Group I. Everything so far has been one sector at a time; here they meet,
and then the result is pushed one stage further — through the **lowering** pass that turns a
polynomial into straight-line C++. By the end of this page you will have seen, in one program, the
entire path a generated kernel takes.

### The diagram

The one-loop quark self-energy: a quark of momentum $p$ emits a gluon, propagates with momentum
$q$, and reabsorbs it. Projected onto the wave-function structure, the Dirac numerator is

$$
T_{\text{num}} \;=\; \mathrm{tr}\!\big[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu\big]\,
P_{\mu\nu}(l), \qquad q = l - p .
$$

Read it as the composition of the two previous steps. The trace is a `DiracNet`
([step-03](step-03.md)) with two open legs $\mu,\nu$ — the two ends of the gluon line. The
projector is an `NNet` ([step-04](step-04.md)) that ties those legs off, carrying the gluon's
transversality and its $1/l^2$. The colour factor $C_F$ ([step-02](step-02.md)) multiplies the whole
thing and is a separate scalar; it does not participate in the index contraction at all.

That factorisation — Dirac chain × Lorentz net × colour number × scalar dressings — is not specific
to this diagram. It is the shape *every* one-loop flow has, and the shape the code generator
assumes.

### The closed form

For $q = l - p$ in the one-angle frame, with $c = \cos\theta$ and $l_1 = |l|$:

$$
T_{\text{num}} \;=\; 4\,p\,\big(-3\,c\,l_1 + p + 2\,c^2 p\big).
$$

Having an independent closed form is the point of the exercise: it pins down the engine's algebra,
including the sign conventions, the $\mathrm{tr}\,\mathbb{1} = 4$, and the Euclidean metric, all at
once.

### Lowering

A contracted diagram is an `MPoly`: a sum of monomials over the frame's symbols and inverse atoms.
Evaluating it as written would mean recomputing shared subexpressions over and over, for every one
of the hundreds of thousands of quadrature points the integrator visits. So the generator does two
build-time passes:

1. **Horner factoring** (`codegen/lower.hpp`) — rewrite the monomial set as a nest of
   `pivot * (…) + (…)` so that shared factors are computed once.
2. **Real value-numbering / CSE** (`codegen/real_cse.hpp`) — accumulate the factored arithmetic
   into a real SSA form that deduplicates every repeated subexpression and folds away trivial
   operations.

The output is a flat `double f(const double*)` with no branches, no divisions, and no temporaries
beyond named scalars. This program runs both passes and prints the result, so you can read the
kernel that a real flow would ship.

## The commented program

`Tutorials/step-05-full-diagram/full_diagram.cpp` does the calculation three times: numerically,
in closed form, and symbolically-then-lowered.

### The frame, as plain numbers

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

For the first pass the momenta are literal numbers, so `nsym = 0` and every component is a constant
polynomial. This is the fastest way to check a value at one kinematic point; there is no symbol
bookkeeping to get wrong.

Note that $q$ is given its own vector id here rather than being written as the combination
$l - p$. Both work. Materialising $q$ is clearer when the components are numbers; the combination
form is what a generated kernel uses, because the frame only carries independent momenta.

### The two networks, contracted together

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: tokens"
:end-before: "@snip end: tokens"
```

This is the first call where both arguments to `numeric_value` are non-empty, and it is worth being
explicit about the order of operations inside. The engine first closes the Dirac chain into
$4\times4$ products and traces it, leaving a tensor with free indices $\mu$ and $\nu$; then it
contracts that tensor against the projector, summing $\mu$ and $\nu$ away; then it returns the
surviving scalar. The $2^{n_p}$ projector expansion and the $(2n-1)!!$ Wick expansion are both
avoided — neither intermediate is ever formed.

### Lowering to a kernel

The first pass used numeric components, so its trace is a constant and would lower to a bare
number. To get a real *kernel* the contraction is redone with the components symbolic:

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: lower"
:end-before: "@snip end: lower"
```

The same `chain` and `lor` objects are reused verbatim — they refer to momentum *ids*, not values,
so they are frame-independent. Only the component table changes. `to_genprog` runs both lowering
passes into a `GenProg` over a shared symbol environment `g`; `emit_env_layout` prints which slot
of the runtime array holds what, and `emit_cpp` prints the function.

## Results

```bash
cmake --build build --target full_diagram && ./build/full_diagram
```

```text
numeric  T_num = 3.527568
closed   T_num = 3.527568   (4 p(-3 c l1 + p + 2 c^2 p), q = l - p)

lowered kernel (symbolic components |p|, l0, l1y):
// fundamental-symbol env layout (fill f[i] per call):
//   f[0] = var(0)
//   f[1] = var(1)
//   f[2] = var(2)
//   f[3] = inv(0)
static inline double T_num([[maybe_unused]] const double *f) {
  const double s0 = -4;
  const double s1 = 4;
  const double s2 = f[0];
  const double s3 = s1*s2;
  const double s4 = f[1];
  const double s5 = s0*s4;
  const double s6 = s3+s5;
  const double s7 = f[2];
  const double s8 = s7*s7;
  const double s9 = s0*s8;
  const double s10 = s4*s6;
  const double s11 = s9+s10;
  const double s12 = -8;
  const double s13 = f[3];
  const double s14 = s11*s13;
  const double s15 = s12+s14;
  const double s16 = s0*s13;
  const double s17 = 8;
  const double s18 = s8*s16;
  const double s19 = s17+s18;
  const double s20 = s2*s19;
  const double s21 = s4*s15;
  const double s22 = s20+s21;
  const double s23 = s2*s22;
  return s23;
}

symbolic T_num = 3.527568   (eval of the lowered polynomial at the same point)
ALL TESTS PASSED
```

### Reading the emitted kernel

This is the artefact the whole library exists to produce, so it repays a careful look.

**The env layout is the calling convention.** `f[0..2]` are the three frame symbols
($|p|$, $l_0$, $l_1$) and `f[3]` is `inv(0)` — the *reciprocal* $1/l^2$. A real kernel's `fill()`
function computes all four from the integrator's arguments once per call. Note that the inverse
atom got a slot of its own: the division happens once, in `fill`, and the trace function itself is
division-free.

**There are 24 operations and no repeats.** Look at `s8 = s7*s7` — that is $l_1^2$, computed once
and used twice (in `s9` and `s18`). That is the CSE pass. Look at the nesting
`s23 = s2*(s20 + s21)` where `s20` and `s21` each contain further products — that is Horner. The
polynomial had a handful of monomials; the emitted form shares every common factor between them.

**Everything is `const double`.** No arrays, no loops, no branches, no function calls. This is what
"straight-line" means, and it is why the same emitted code compiles unchanged for a GPU
([step-21](step-21.md)): there is nothing in it a device compiler could object to.

**The constants are folded.** `s0 = -4` and `s17 = 8` appear as literals because the coefficients
were exact rationals in the polynomial. Nothing is computed at runtime that could be computed at
build time.

### And the value is right, three ways

The numeric contraction, the closed form, and an evaluation of the lowered polynomial all give
3.527568. The third check matters more than it looks: it is the guarantee that lowering is
*value-preserving*. A Horner or CSE bug would show up here and nowhere else, since the first two
paths never touch `to_genprog`.

## Possibilities for extensions

1. **Add the colour factor.** Multiply by the $C_F = 4/3$ from [step-02](step-02.md) and confirm
   the emitted kernel changes by exactly one constant. In a real flow the colour factor is folded
   into the diagram's scalar coefficient, never into the trace — see if you can see why from the
   emitted code.

2. **The second diagram.** The real quark self-energy has *two* diagrams: the regulator insertion
   sits on the gluon line in one and on the quark line in the other, which amounts to $q = l - p$
   versus $q = l + p$. Redo the contraction with the sign flipped and compare the polynomials. They
   should share most of their structure — which is exactly the redundancy the generator's
   cross-diagram deduplication exploits ([step-20](step-20.md)).

3. **Watch CSE earn its keep.** Comment out the projector's atom (use a plain metric instead) and
   re-emit. Then put it back but add a second projector. Count the emitted lines each time. The
   growth is sublinear in the monomial count, and that gap is the whole value of the pass.

4. **Change the frame and re-emit.** Give $l$ a third nonzero component. The value at a fixed point
   must not change; the emitted kernel will grow. This is [step-04](step-04.md) extension 6 seen
   from the other end — you are now looking at what the extra components *cost*.

5. **Break it deliberately — wrong atom value.** Pass `l2` instead of `1.0/l2` to `eval`. You will
   get a wrong number with no complaint whatsoever. Sit with that for a moment: the atom convention
   is a contract between you and the engine, and it is not checkable. This is precisely why real
   kernels are generated rather than hand-written, and why [step-22](step-22.md) is a whole page
   about validation.

6. **Read the lowering source.** With the emitted output in front of you,
   [CSE and Horner lowering](../internals/cse-and-lowering.md) will make considerably more sense
   than it would have before.

## The plain program

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
```

You have now done by hand everything the code generator does. [step-06](step-06.md) shows the
front-end doing all of it from a one-line description — and the kernel it emits will look extremely
familiar.
