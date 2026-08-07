# step-04: Lorentz networks and inverse atoms

*Builds on: [step-01](step-01.md) · Built on by: [step-05](step-05.md), [step-16](step-16.md) ·
Tags: `lorentz`, `projector` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

[step-01](step-01.md) contracted a projector with numeric components. This step does it with the
components **symbolic**, which changes the character of the result completely: instead of a number
you get a polynomial, and that polynomial is what a generated kernel evaluates. Two ideas appear
for the first time — the *frame*, and the *inverse atom*.

### The frame

A loop integrand depends on the external momenta and the loop momentum, but not on all of their
components independently: rotational invariance means only magnitudes and relative angles matter.
A **frame** is the choice of concrete components that encodes this. For a one-loop propagator
diagram — one external momentum $p$, one loop momentum $l$ — the natural choice is

$$
p = (p, 0, 0, 0), \qquad l = (l_0, l_1, 0, 0),
$$

i.e. put $p$ along an axis and let $l$ live in the plane it spans with that axis. Two components of
$l$ are identically zero and can be dropped; the integrand is a function of exactly three scalars.

This is worth stating carefully because it is the single biggest lever on kernel size. The engine
contracts *over the frame's components*. Choosing a frame with two zero components means those
entire index-sum branches vanish before any polynomial is built — the saving is not a
simplification afterwards, it is work never done. [step-07](step-07.md) is entirely about frames
and the five builders the front-end provides; here we hand-write one.

### The inverse atom

A transverse projector

$$
P(l)_{\mu\nu} = \delta_{\mu\nu} - \frac{l_\mu l_\nu}{l^2}
$$

carries a denominator. In a symbolic frame $l^2 = l_0^2 + l_1^2$, a polynomial — and dividing one
polynomial by another does not in general give a polynomial. So the engine cannot simply "do the
division"; it must carry $1/l^2$ as an opaque quantity.

It does this by giving each distinct denominator an **atom id**, and letting monomials carry
powers of atoms alongside powers of symbols. A monomial is then

$$
c \cdot \prod_i x_i^{e_i} \cdot \prod_j \mathrm{atom}_j^{a_j},
$$

and multiplication just adds exponents. Crucially the engine is *also* told what each atom is the
reciprocal of, which lets it cancel: when a monomial acquires both an $l^2$ from a numerator and an
$\mathrm{atom}_0$, the two annihilate instead of both being carried. Without that cancellation the
polynomials would grow without bound through a long chain of propagators.

```{admonition} Why the caller numbers the atoms
:class: note
Atom ids are supplied by you, not invented by the engine, and that is deliberate: two projectors on
the *same* momentum must share an atom id so their denominators cancel against the same numerators.
In generated kernels this bookkeeping is done by the front-end, which allocates one atom per
distinct inverse propagator across the whole diagram. When hand-building, it is your job — and
giving two projectors on the same momentum *different* ids is a correctness bug that no test will
catch for you, only a size regression.
```

## The commented program

`Tutorials/step-04-lorentz-networks/lorentz_networks.cpp` computes $p\cdot P(l)\cdot p$ in the
one-angle frame.

### Writing the frame down

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

Three symbols for three nonzero components. Everything else is `env.zero()` — and note that this is
a *structural* zero, not the number 0.0: a zero `MPoly` has no monomials, so any product it enters
is dropped immediately rather than being carried as a term with coefficient zero.

Compare this with [step-03](step-03.md), which used eight independent symbols for two momenta. Both
are legitimate; they answer different questions. Eight symbols verify an identity *in general*;
three symbols compute the thing a kernel actually needs. Real generation always uses the frame.

### Declaring the atom

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
:start-after: "@snip begin: atom"
:end-before: "@snip end: atom"
```

`atomDen[0]` is $l^2$ as a polynomial, built from the same `comp` table the network will use — so
it cannot disagree with the frame. This is the reciprocal-of information from the introduction:
atom 0 means $1/l^2$, and the engine now knows what to cancel it against.

### The network, and the contraction

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
:start-after: "@snip begin: net"
:end-before: "@snip end: net"
```

`nprojT(mu, nu, {{1.0, 1}}, 0)` is the projector: legs `mu` and `nu`, built on momentum id 1
(the loop momentum $l$), with denominator atom 0. Both labels appear twice across the three
factors, so both are summed and the network closes.

The Dirac argument is the empty `DiracNet{}` — this network is pure Lorentz. Every combination is
allowed here: pure Lorentz (this step), pure Dirac ([step-03](step-03.md) part a), or both together
([step-05](step-05.md)).

## Results

```bash
cmake --build build --target lorentz_networks && ./build/lorentz_networks
```

```text
contracted monomials = 2   (p.P.p = sp(p,p) - sp(p,l)^2 / l^2)
p.P(l).p             = 1.11905   (= p^2 (1 - cos^2) = 1.11905)
p.P(l).p / p^2       = 0.662162   (= 1 - cos^2 theta = 0.662162)
ALL TESTS PASSED
```

**Two monomials.** That is the headline. The contraction ran over $\mu,\nu \in \{0,1,2,3\}$ —
sixteen index combinations — and what came back is

$$
p_0^2 \;-\; p_0^2\,l_0^2 \cdot \mathrm{atom}_0 ,
$$

two terms. The frame did the work: because $p$ has only component 0, every term involving $p_1$,
$p_2$, $p_3$ was structurally absent. This is the compactness that makes generated kernels small.

**The physics.** $p\cdot P(l)\cdot p / p^2 = 1 - \cos^2\theta$, where $\theta$ is the angle between
$p$ and $l$. This is the angular weight of a transverse gluon exchange — the factor that appears in
the ghost and gluon loops of every Yang–Mills propagator flow. In a real kernel it is multiplied by
dressings and a regulator and handed to a quadrature over $\cos\theta$.

**The asymmetry in the atom, again.** `atomDen` got $l^2$; `eval` got `1.0/l2v`. The engine wants
the denominator to reason about; the evaluator wants the reciprocal precomputed. A generated kernel
computes each `1/l²` once per call and then never divides again — which is why the emitted
arithmetic is division-free and GPU-friendly.

## Possibilities for extensions

1. **Watch the cancellation happen.** Contract $l\cdot P(l)\cdot l$, which is analytically zero.
   Check that `poly.size() == 0` — the two monomials cancel *during* contraction because the
   numerator $l^2$ met atom 0. Now deliberately give the projector a fresh atom id 1 with the same
   denominator, and observe that you get two monomials that only cancel numerically at `eval` time.
   That is the bug the "share atom ids" rule prevents.

2. **Two projectors.** Add a second `nprojT` on the same momentum sharing a middle label, i.e.
   $p\cdot P(l)P(l)\cdot p$. Idempotence says the answer is unchanged. Check both the value *and*
   the monomial count — if the count doubled, your atom ids are wrong.

3. **The longitudinal complement.** Compute $p\cdot P^L(l)\cdot p$ with `nprojL` and confirm it
   equals $p^2\cos^2\theta$, and that the two add to $p^2$.

4. **Change the frame.** Put $p$ along axis 1 instead of axis 0 and re-run. The value must not
   change (it is a scalar), but look at the polynomial — the *monomials* may differ. Then try
   giving $l$ three nonzero components and watch the monomial count grow. This is the frame lever
   from the introduction, made visible.

5. **Break it deliberately — a free index.** Drop the second `nvec`, leaving `nu` open. The network
   no longer closes. Read the error and remember what it looks like; an unclosed index in a
   hand-built net is one of the two mistakes ([step-01](step-01.md) extension 3 is the other) that
   account for most first-day confusion.

6. **A frame with no zero components.** Use four symbols for $l$ and four for $p$ and redo the
   contraction, then compare the monomial count to the three-symbol version. The ratio is the price
   of not choosing a frame — and for a four-point vertex it is the difference between a kernel that
   generates and one that does not ([step-20](step-20.md)).

## The plain program

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
```

Next, [step-05](step-05.md): Dirac and Lorentz together — a whole diagram — and then the lowering
pass that turns the resulting polynomial into C++.
