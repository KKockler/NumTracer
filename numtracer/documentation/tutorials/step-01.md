# step-01: Hello, tensor network

*Builds on: nothing · Built on by: [step-02](step-02.md), [step-03](step-03.md),
[step-04](step-04.md) · Tags: `network`, `projector` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

NumTracer exists to answer one question over and over: *given a product of tensors with shared
indices, what single number does it collapse to?* Everything else in the library — the code
generator, the Dirac algebra, the SU($N$) tables, the CSE pass — is scaffolding around that one
operation. So the first program does exactly it, with no physics in sight.

### The one rule

A **network** is a product of tensor factors, each carrying integer *index labels*. Two index slots
are summed over exactly when they carry **the same label** — the Einstein convention, made literal.
A label that appears once stays free; a network with no free labels left is a scalar.

That rule is the entire mental model, and it is worth being precise about what it does *not* say.
It says nothing about "upper" and "lower" indices — NumTracer is Euclidean, so raising and lowering
is the identity and there is no distinction to track. It says nothing about the *extent* of an
index: a Lorentz label runs over 4 values, a fundamental SU(3) label over 3, a spinor label over 4,
and the contraction machinery never needs to know which is which, because it only ever compares
labels for equality. This is why a single network can mix all three sectors without them
interfering (see [Scope & conventions](../getting_started/scope-and-conventions.md)).

### The two networks

We contract two networks and check each against its closed form. First, the plainest possible
contraction — two vectors tied through a metric:

$$
a\cdot b \;=\; a_\mu\,\delta^{\mu\nu}\,b_\nu .
$$

Written as a network that is three factors — $a$ on label $\mu$, the metric on $\mu$ and $\nu$,
$b$ on $\nu$ — and both labels are used twice, so both are summed away.

Second, something with structure: a **transverse projector**,

$$
P(k)_{\mu\nu} \;=\; \delta_{\mu\nu} - \frac{k_\mu k_\nu}{k^2},
\qquad\text{so}\qquad
a\cdot P(k)\cdot a \;=\; a^2 - \frac{(a\cdot k)^2}{k^2}.
$$

$P(k)$ is the object that projects out the component of a vector along $k$; in a gauge theory it is
what sits on an internal gluon line. It matters here for a structural reason: it carries a
**denominator**, $1/k^2$. Contraction cannot evaluate that — $k^2$ is generally a runtime quantity —
so the engine must carry it through symbolically. The way it does so (an *inverse atom*) is the
subject of [step-04](step-04.md); this program just meets the idea.

```{admonition} Why not simply build the tensors and sum?
:class: note
You could materialise $P(k)$ as sixteen numbers and loop. That works here and does not work later:
a Dirac trace of $n$ gammas has $(2n-1)!!$ index pairings, and a Lorentz trace of $n_p$ transverse
projectors has $2^{n_p}$ terms if you expand the definition. The engine's contraction never forms
either intermediate. Step 1 uses the small case to introduce the vocabulary; steps 2–4 show the
same vocabulary doing the thing you could not do by hand.
```

## The commented program

The whole program is `Tutorials/step-01-hello-network/hello_network.cpp`. We take it in three
pieces.

### Labels and vectors

Index labels are just integers, but writing raw numbers invites accidental collisions — two factors
that were never meant to contract quietly sharing a label is a silent wrong answer, not a
compile error. The idiom used throughout the tutorials is a single unscoped `enum`: every name gets
a distinct auto-numbered value, so factors contract only where you deliberately reuse a name.

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:lines: 24-33
```

Note that the *vector* ids (`a`, `b`, `k` as momenta) are a separate namespace from the *index*
labels. A vector id indexes the component table below; an index label names a slot to be summed.
Confusing the two is the most common beginner error.

### The frame

Before anything can be contracted, every vector needs components. The table that supplies them is
called the **frame**:

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

`comp[vid][c]` is component `c` (0–3) of vector `vid`. The entries are `MPoly` — multivariate
polynomials — not `double`, because in a real kernel the components are *symbolic*: functions of
the loop magnitude and angles that the integrator supplies at runtime. Here they happen to be
literal numbers, which is the degenerate case `nsym = 0`: an empty symbol space, every component a
constant polynomial.

`LorentzEnv env(nsym)` is the factory every polynomial is minted through. It exists to bind `nsym`
**once**. Two polynomials can only be added or multiplied if they agree on how many symbols the
world has, and mixing two different `nsym` values is undefined behaviour rather than a checked
error — so the library makes the size unforgeable by hanging every constructor off an env that
already knows it. [step-03](step-03.md) returns to this.

### Network one: a metric

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: dot"
:end-before: "@snip end: dot"
```

An `NNet` is a **sum of terms**; an `NTerm` is a coefficient times a product of factors. Here there
is one term, coefficient 1, three factors. `nvec(mu, {{1.0, 0}})` is "vector id 0, coefficient 1, on
index `mu`" — the `{{coeff, vid}}` list is a *linear combination*, which is what lets a propagator
momentum $l - p$ be written `{{1.0, l}, {-1.0, p}}` without introducing a new vector.

`numeric_value` is the engine's one contraction entry point. Its first argument is the Dirac chain,
empty here because this network is pure Lorentz. It returns an `MPoly`; `eval` then substitutes
numbers — symbol values (none) and inverse-atom values (none).

### Network two: a projector

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: proj"
:end-before: "@snip end: proj"
```

Two things are new. `nprojT(mu, nu, {{1.0, 2}}, 0)` takes a fourth argument, `0`: the **id of the
inverse atom** that holds this projector's $1/k^2$. Atoms are numbered by the caller, and the
engine needs to know what each one is the reciprocal *of* — hence `k2`, computed just above and
handed to `numeric_value` as the atom-denominator table. Knowing the denominator is what lets the
engine cancel a $k^2$ against a $1/k^2$ inside a monomial instead of carrying both.

Second, `eval` now receives `{1.0 / k2v}`: one value per atom. Note the asymmetry — the *engine*
is told the denominator ($k^2$), the *evaluator* is told the reciprocal ($1/k^2$). That is not an
inconsistency: the engine needs $k^2$ to do algebra with, while the evaluator wants the reciprocal
precomputed so the generated kernel never divides.

## Results

```bash
cmake --build build --target hello_network && ./build/hello_network
```

```text
a . b        = 0.26   (expect 0.26)
a . P(k) . a = 0.403649   (expect a^2 - (a.k)^2/k^2 = 0.403649)
ALL TESTS PASSED
```

Both agree with the closed form to machine precision, which is the point: the engine was told
nothing about what a metric or a projector *means*. It was told that $\delta_{\mu\nu}$ has entries
$1$ on the diagonal and that $P$ has entries $\delta_{\mu\nu} - k_\mu k_\nu \cdot (\text{atom }0)$,
and it summed the shared labels. The identity $a\cdot P(k)\cdot a = a^2 - (a\cdot k)^2/k^2$ is not
knowledge the engine has — it is what falls out.

Notice also what the second result *is*. With $k$ fixed and $a$ varying, $a\cdot P(k)\cdot a / a^2$
is $1 - \cos^2\theta$, the angular weight that appears in every one-loop self-energy. You have
already computed a piece of real physics; [step-04](step-04.md) does it with the angle symbolic.

## Possibilities for extensions

1. **Check that the projector projects.** $P(k)$ should annihilate $k$: contract
   `nvec(mu, {{1.0, 2}}) · nprojT(mu, nu, {{1.0,2}}, 0) · nvec(nu, {{1.0, 2}})` and confirm you get
   0 to machine precision. Then check idempotence, $P\cdot P = P$, by tying two projectors through
   a shared middle label and comparing to a single one.

2. **A momentum that is a combination.** Replace the second `a` in network two with $a - b$,
   i.e. `nvec(nu, {{1.0, 0}, {-1.0, 1}})`. Predict the closed form first, then check it. This is
   the mechanism every propagator momentum uses.

3. **Break it deliberately — reuse a label.** Change the metric to `nmet(mu, mu)`. You will not get
   an error; you will get $\mathrm{tr}\,\delta = 4$ times something, because a label used twice
   *within one factor* is a self-contraction and is perfectly legal. Now instead give the two
   `nvec` factors the *same* label as each other and no metric at all, and reason about what you
   get. The lesson is that the engine cannot tell an intended contraction from a typo — labels are
   the whole contract, which is why the `enum` idiom is not decoration.

4. **Break it deliberately — forget the atom.** Pass `{}` instead of `{k2}` as the atom-denominator
   table, or `{}` instead of `{1.0/k2v}` to `eval`. Read what happens carefully: one of these is a
   loud failure and one is a quiet wrong number. Knowing which is which will save you an afternoon
   later.

5. **Add a longitudinal projector.** `nprojL` is $k_\mu k_\nu / k^2$. Verify $P^T + P^L = \delta$
   by contracting all three against the same pair of vectors and comparing sums.

## The plain program

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
```

Next, [step-02](step-02.md): the same idea in the SU($N$) sector, where the network folds to a
number without the tensor ever existing.
