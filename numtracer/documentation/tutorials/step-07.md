# step-07: Frames and kinematics

*Builds on: [step-06](step-06.md) · Built on by: [step-10](step-10.md), [step-16](step-16.md) ·
Tags: `frame`, `kinematics` · **Tier B** (a Wolfram kernel and FunKit)*

## Introduction

The frame is the least glamorous part of a flow and the part most likely to cost you a week. It is
where physics kinematics meets the engine's component arithmetic, and a frame that is subtly not
the kinematics you meant produces a kernel that compiles, runs, and is wrong.

### A frame is a component table

That is the whole definition. A frame is an Association mapping each momentum tag to its four
components, written in terms of the kernel's runtime scalars. `propFrame` produces:

```text
q1 = {p, 0, 0, 0}
ql = {cos1*l1, Sqrt[1 - cos1^2]*l1, 0, 0}
```

That is it. The external momentum sits along axis 0 with magnitude `p`; the loop momentum has
magnitude `l1` at angle `cos1` to it, in the plane those two span. Two components of `ql` are
identically zero and every index sum over them vanishes structurally — which, as
[step-04](step-04.md) showed, is where kernel compactness comes from.

**Printing a frame is the single most useful debugging habit in this whole toolchain.** It is data;
look at it.

### Why the choice is constrained

A loop integrand is a scalar, so it can only depend on rotational invariants: magnitudes and
relative angles. With $n$ external momenta plus one loop momentum in 4 dimensions, the number of
independent invariants is fixed, and a frame is a *parametrisation* of them. Two requirements pull
against each other:

* **Completeness.** The frame must span every configuration the integrand is evaluated at. A frame
  that accidentally imposes a relation between momenta silently restricts the calculation.
* **Compactness.** Every nonzero component costs contraction work and kernel size. You want as many
  structural zeros as the physics allows.

The builders in `Frames.m` are the standard resolutions of that tension for the kinematics fRG
flows actually use.

### The five builders

| Builder | Kinematics | Loop angles | Typical use |
|---|---|---|---|
| `propFrame[p, l1, cos1, q1, ql]` | 1 external | 1 | propagator flows: $Z_A$, $Z_c$, $Z_q$ |
| `sp3Frame[p, l1, cos1, cos2, q1, q2, q3, ql]` | 3 externals at the **symmetric point** (120° apart, equal magnitudes) | 2 | 3-point vertices: $Z_{A^3}$, $Z_{Ac\bar c}$ |
| `gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, …]` | 3 externals, **general** (independent magnitudes and opening angle) | 2 | 3-point vertices where the symmetric point is degenerate |
| `sp4Frame[p, l1, cos1, cos2, phi, q1…q4, ql]` | 4 externals, regular tetrahedron | 3 | 4-point vertices: $Z_{A^4}$, four-quark |
| `propFrameFT[p0, p, l0, l1, cos1, pSym, lSym]` | finite $T$: independent temporal components | 1 spatial | thermal flows ([step-16](step-16.md)) |

### The trap: the symmetric point can be degenerate

`sp3Frame` is compact and standard, and it is what you should reach for first. But it is a *slice*
of the configuration space, and some tensor bases are **degenerate on that slice** — their Gram
matrix has vanishing determinant there, so the projectors that invert it do not exist. The full
`AqbqDirect` quark–gluon vertex basis is one such case.

When that happens you have not found a bug; you have found that your basis needs more kinematic
freedom than the symmetric point provides. `gen3Frame` supplies it, at the cost of three runtime
arguments instead of one and a larger kernel.

This step's self-check is precisely the relationship between the two, because "the general frame
reduces to the special one where both are valid" is what licenses trusting the general one where
only it is.

## The commented program

`Tutorials/step-07-frames/frames.wls`.

### Looking at the frames

```{literalinclude} ../../../Tutorials/step-07-frames/frames.wls
:language: mathematica
:start-after: "@snip begin: show"
:end-before: "@snip end: show"
```

### One network, two frames

```{literalinclude} ../../../Tutorials/step-07-frames/frames.wls
:language: mathematica
:start-after: "@snip begin: nets"
:end-before: "@snip end: nets"
```

The network mentions only momentum *tags* and index labels — it contains no kinematics at all. That
is the separation the DSL is built around: the network is the algebra, the frame is the kinematics,
and you can swap either without touching the other. Everything below reuses `net` verbatim.

```{admonition} Those parentheses are load-bearing
:class: warning
In a `.wls` script a line that is already a complete expression **ends there**, so

```mathematica
net = termA
    + termB;
```

assigns only `termA` and silently discards the rest — producing a kernel that compiles, runs, and
computes the wrong network. This script hit exactly that while being written. The guard is the
`expected 2 diagrams` assertion in the emit helper below; assert the count you expect, not merely a
nonzero one. [step-09](step-09.md) discusses it further.
```

Note the `"Args"` lists differ. `sp3Frame` needs four runtime scalars, `gen3Frame` six. The `"Args"`
list is exactly the generated C++ signature, in order.

## Results

```bash
cmake --build build --target frames && ./build/frames
```

The script first prints the three frames:

```text
--- propFrame[p, l1, cos1, q1, ql]  (1 external, 1 angle) ---
  q1 = {p, 0, 0, 0}
  ql = {cos1*l1, Sqrt[1 - cos1^2]*l1, 0, 0}
--- sp3Frame[p, l1, cos1, cos2, q1, q2, q3, ql]  (3 externals at the symmetric point, 2 angles) ---
  q1 = {p, 0, 0, 0}
  q2 = {-1/2*p, (Sqrt[3]*p)/2, 0, 0}
  q3 = {-1/2*p, -1/2*(Sqrt[3]*p), 0, 0}
  ql = {cos1*l1, Sqrt[1 - cos1^2]*cos2*l1, Sqrt[(1 - cos1^2)*(1 - cos2^2)]*l1, 0}
--- gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, q1, q2, q3, ql]  (3 externals, general) ---
  q1 = {p1m, 0, 0, 0}
  q2 = {cosP*p2m, Sqrt[1 - cosP^2]*p2m, 0, 0}
  q3 = {-p1m - cosP*p2m, -(Sqrt[1 - cosP^2]*p2m), 0, 0}
  ql = {cos1*l1, Sqrt[1 - cos1^2]*cos2*l1, Sqrt[(1 - cos1^2)*(1 - cos2^2)]*l1, 0}
```

Read these three tables carefully; they contain most of what there is to know.

* **`sp3Frame` hard-codes the symmetric point.** $q_2$ and $q_3$ have components
  $(-\tfrac12 p, \pm\tfrac{\sqrt3}{2}p)$ — that is 120° apart, all three of equal magnitude $p$.
  The whole configuration is one number.
* **`gen3Frame` does not.** $q_1$ has magnitude `p1m`, $q_2$ has magnitude `p2m` at angle `cosP`,
  and $q_3 = -(q_1 + q_2)$ by momentum conservation. Substituting `p1m = p2m = p`,
  `cosP = -1/2` turns the third table into the second — $\cos 120° = -1/2$ — which is the identity
  the checker verifies.
* **Both 3-point frames give the loop momentum three nonzero components**, because two angles are
  needed to orient it relative to a plane. `propFrame` needs only two components. That extra
  component is most of why a 3-point kernel is an order of magnitude bigger than a propagator one.
* **The fourth component is always zero.** With three external momenta in a plane plus a loop
  momentum, nothing reaches out of the 3-space they span. `sp4Frame` is where the fourth component
  finally switches on — and where the third loop angle `phi` appears.

Then the checker:

```text
step-07: sp3_kernel diagrams = 2
step-07: gen3_kernel diagrams = 2

sp3Frame  vs  gen3Frame at (p1m = p2m = p, cosP = -1/2)
  1296 kinematic points, worst |difference| = 3.553e-15   (typical |value| ~ 10)
  points disagreeing beyond 1e-11 relative : 0
  off the symmetric point they differ      : yes (-1.40568 vs -0.954884)
ALL TESTS PASSED
```

The second check is not decoration. A test that only ever evaluates the general kernel *at* the
symmetric point would pass just as happily if `gen3Frame` ignored its extra arguments entirely. The
off-point comparison is what makes the on-point agreement mean something — a habit worth copying
into your own gates ([step-22](step-22.md)).

## Two things you will meet in a real flow

### `AngleDefs`: keeping the scalar and tensor halves consistent

A multi-angle flow has two halves that must agree about kinematics. The *tensor* half gets its
angles from the frame. The *scalar* half — propagator denominators, dressings evaluated at
$\sqrt{(l-p)^2}$ — is reduced separately in Mathematica, and it naturally produces symbolic
loop–external cosines like `cosl1p1`, `cosl1p2`, `cosl1p3`.

`"AngleDefs"` is the map that ties them together. In `QCD_Nf2/baseline/QCD.wl`:

```mathematica
spAngles3 = {
  cosl1p1 -> cos1,
  cosl1p2 -> (-cos1 + Sqrt[3 - 3 cos1^2] cos2)/2,
  cosl1p3 -> (-cos1 - Sqrt[3 - 3 cos1^2] cos2)/2};
```

Those expressions are exactly $\hat l \cdot \hat q_i$ computed from the `sp3Frame` table above — you
can derive them from the printed components in a line each. **They must match the frame**, and
nothing checks that they do: an `AngleDefs` written for `sp3Frame` and used with `gen3Frame` gives
a kernel whose tensor and scalar halves disagree about where the momenta are.
[step-14](step-14.md) uses them for real.

### Hand-rolled frames

Nothing requires you to use a builder. A frame is an Association, and you can write one. The Yukawa
flow in `QCD_Nf2/with_mesons` does, because it needs the $(S_0, S_1, S_\varphi)$ shape coordinates
that its interpolation grid is defined on, and its author left a warning worth reproducing:

```{admonition} From QCD_Nf2/with_mesons/QCD.wl — why gen3Frame could not be used
:class: warning
> The earlier scaffolding here built the frame from `NumTracer`gen3Frame` and then mapped its
> scalars through `TB3PToS0S1SPhi`. That is NOT consistent: `gen3Frame` spans the plane with `p1`
> along `e0` and `p2` along `+e1`, whereas `TB3PToS0S1SPhi` spans it with `Q = p1` along `e0` and
> `k = (p3-p2)/2` along `+e1` — and `p2 = -k - Q/2`, so the two differ by `e1 -> -e1`.

Two parametrisations of "three momenta in a plane" that are individually correct can still differ
by a reflection, and composing them silently gives the wrong kinematics. If you build a frame by
hand, print it and check a scalar product you know.
```

The same file carries a second warning, about applying a momentum transformation to a whole trace:

```{admonition} From QCD_Nf2/with_mesons/QCD.wl — do not transform a trace
:class: warning
> `transf3PTo` must NOT be applied to a whole trace. Its first act is a GLOBAL momentum
> substitution `p1 -> Q`, … Momenta sitting in TENSOR slots (`pdash[p1,...]`, `ntVec[...]`) are
> rewritten by the first step and then left as the Module-local `Q`, `k` … They survive all the way
> into the generated kernel as free variables named `"Q$21048"`, `"k$21048"` and the generator fails
> to compile.

A `$`-suffixed symbol in a generated file always means a Module-local leaked into emitted code.
Scalar reductions belong on the scalar half only; the tensor half gets its kinematics from the
frame.
```

## Under the hood: the unit-loop fast path

Most vacuum frames have a special structure: every component of the loop momentum is `l1` times
something that depends only on angles. The engine detects this and factors `l1` out of the
contraction entirely, working with a *unit* loop momentum and restoring the power at the end. That
is a large saving, and it is why the vacuum kernels are as small as they are.

A frame that breaks the pattern — `propFrameFT`, whose temporal loop component `l0` is independent
of `l1` ([step-16](step-16.md)) — automatically takes the general path instead, which mints extra
sine/cosine/radical symbols (`polyFrameSpec`). Nothing is wrong when this happens; the kernel is
simply bigger. The environment variable `NT_NO_UNIT_GROUPS` forces the general path, which is
useful exactly once: to check that a suspicious unit-loop kernel agrees with the general one.

## Possibilities for extensions

1. **Derive `AngleDefs` yourself.** From the printed `sp3Frame` table, compute $\hat l\cdot\hat q_2$
   by hand and confirm you get $(-\cos_1 + \sqrt{3-3\cos_1^2}\,\cos_2)/2$. Then do the same for
   `gen3Frame` and note that the answer now involves `cosP` — which is why an `AngleDefs` cannot be
   shared between frames.

2. **Add `sp4Frame`.** Generate the same network on the four-point frame (you will need a fourth
   momentum tag and the `phi` argument). Compare `nenv` and the emitted line count against the
   3-point kernels. This is the growth curve that makes [step-20](step-20.md) necessary.

3. **Verify momentum conservation.** In `gen3Frame`, check that $q_1 + q_2 + q_3 = 0$ componentwise
   from the printed table. Then generate a network containing `ntSP[q3, q3]` and confirm the kernel
   agrees with $(q_1+q_2)^2$ computed from the scalars.

4. **Break it deliberately — mismatched `AngleDefs`.** Generate a network with a scalar coefficient
   containing `cosl1p2`, and supply an `AngleDefs` written for the *wrong* frame. The kernel will
   compile and run. Compare it against the correct one and see how large the error is — this is the
   failure mode that has no error message at all.

5. **Feel the degeneracy.** Take the `gen3Frame` kernel and evaluate it approaching the symmetric
   point along different directions in `(p1m, p2m, cosP)`. For this simple network everything is
   smooth. Then read the `AqbqDirect` note in [step-11](step-11.md) for the case where it is not.

## The plain program

```{literalinclude} ../../../Tutorials/step-07-frames/frames.wls
:language: mathematica
```

```{literalinclude} ../../../Tutorials/step-07-frames/check_frames.cpp
:language: cpp
```

Next, [step-08](step-08.md): the complete DSL head vocabulary — every tensor the engine knows,
each with a network small enough to check by hand.
