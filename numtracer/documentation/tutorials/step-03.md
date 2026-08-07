# step-03: Dirac traces — matrices, then tokens

*Builds on: [step-01](step-01.md) · Built on by: [step-05](step-05.md), [step-18](step-18.md) ·
Tags: `dirac`, `trace` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

The Dirac sector is where the naive approach stops being merely wasteful and becomes impossible, so
it is worth being precise about what a gamma trace *is* before writing any code.

### What a trace actually is

A quark propagating through a diagram carries a spinor index. Every vertex it passes inserts a
$4\times4$ matrix on that index; when the quark line closes into a loop, the first and last spinor
indices are identified and summed. That sum-over-the-diagonal is the trace. So

$$
\mathrm{tr}\big[\slashed p\,\slashed q\big]
$$

is *literally* the trace of a product of two $4\times4$ matrices, where a **slashed momentum**
$\slashed p = p^\mu\gamma_\mu$ is the matrix you get by weighting the four gammas with the four
components of $p$.

The textbook way to evaluate such a trace is Wick's theorem: sum over all pairings of the gamma
indices, with a sign per crossing. For $2n$ gammas that is $(2n-1)!!$ terms — 3 for four gammas,
10395 for fourteen, and a four-point vertex flow reaches well past that. Worse, the intermediate is
a sum of scalar-product monomials that must then be collected.

NumTracer does not use Wick's theorem. It **multiplies the matrices**. The gammas are typed out as
`constexpr` tables in `dirac/dirac_data.hpp` (Hermitian, chiral/Weyl basis, Euclidean, so
$\{\gamma^\mu,\gamma^\nu\} = 2\delta^{\mu\nu}$), each with only four nonzero entries out of
sixteen. Multiplying a chain is $O(n)$ matrix products of a fixed $4\times4$ size, and the Clifford
algebra is *implied* rather than applied — every identity you know falls out of the arithmetic.

```{admonition} The entries are polynomials, not numbers
:class: important
This is the trick that makes the matrix-product approach work symbolically. The engine does not
multiply matrices of `double`; it multiplies matrices of `MPoly` — multivariate polynomials in the
frame's symbols. So $\slashed p$ is a $4\times4$ array whose entries are things like
"$p_0 + i p_3$", and the product of two such matrices has entries that are polynomials of degree 2.
The trace is then one polynomial, already collected, with no separate simplification pass.
```

This step shows the same trace twice: first with the matrices explicit, so you can see there is no
magic; then through the **token API**, which is how every generated kernel expresses it.

## The commented program, part a — as a matrix product

`Tutorials/step-03-dirac-traces/trace_raw.cpp`.

### The symbol space

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.cpp
:language: cpp
:start-after: "@snip begin: symbols"
:end-before: "@snip end: symbols"
```

Here `nsym = 8` and the components are *bare symbols*: $p_\mu = x_\mu$, $q_\mu = x_{4+\mu}$. This
is the fully general case — no frame has been chosen, the eight components are independent — and it
is the right setting to check an algebraic identity, because if the polynomial identity holds for
independent symbols it holds in every frame.

```{admonition} `nsym` is a dimension you fix up front, not a running count
:class: important
The `nsym` handed to `LorentzEnv env(nsym)` is **not** "how many variables I have created so far".
It is the fixed dimension of the variable space: internally every monomial stores an exponent
vector of length `nsym`, one slot per symbol. You choose it once, then *address* symbols by index —
`env.var(0)` is $x_0$, `env.var(7)` is $x_7$. Symbol `i` exists for any `0 ≤ i < nsym` whether or
not you use it, and its value is supplied later, positionally, in the `eval` array.

**Every `MPoly` that is combined must share the same `nsym`.** `+` and `*` walk both operands'
exponent vectors slot-for-slot; mixing two different sizes reads out of bounds — undefined
behaviour, not a checked error — and like terms silently fail to combine. Building everything from
**one** env guarantees this. When momenta are plain numbers rather than symbols (as in
[step-05](step-05.md)) you make a `LorentzEnv env(0)` and use `env.constant`.
```

### Slash, multiply, trace

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.cpp
:language: cpp
:start-after: "@snip begin: slash"
:end-before: "@snip end: slash"
```

Three lines, and every one is ordinary linear algebra. `env.slashC(comp)` builds
$\sum_\mu \mathtt{comp}[\mu]\gamma^\mu$ as a `Mat4` of `MPoly`; `matmul` is the $4\times4$ product;
`mtrace` sums the diagonal. There is no Dirac-algebra code path here at all — `matmul` does not
know it is multiplying gammas.

## Results, part a

```bash
cmake --build build --target trace_raw && ./build/trace_raw
```

```text
tr(p/ q/)  = 1.04 + 0i   (4 monomials in the polynomial)
4 (p.q)    = 1.04
ALL TESTS PASSED
```

Two things to notice.

**The identity was never coded.** NumTracer contains no rule saying
$\mathrm{tr}(\slashed p\slashed q) = 4\,p\cdot q$. It multiplied two matrices of polynomials and
traced. The result agrees because the tabulated gammas satisfy the Clifford algebra — which is
itself checked, in `tests/`, against a runtime oracle, so a typo in the tables fails a test rather
than silently producing wrong physics.

**"4 monomials" is the whole answer.** The polynomial is $4\sum_\mu p_\mu q_\mu$ — four terms, one
per $\mu$, each with coefficient 4. That *is* $4\,p\cdot q$, written in components. And the
imaginary part is exactly zero, not $10^{-17}$: the individual matrix entries are complex (the
gammas have $\pm i$ in them) and the imaginary parts cancel term by term in exact arithmetic.

## The commented program, part b — via the token API

In practice you never build the matrices. You **describe** the closed chain as a list of tokens and
let the engine contract it — which is exactly the form the code generator emits.

```{admonition} What `numeric_value` does
:class: note
`env.numeric_value(dirac, lorentz, comp, atomDen)` is the engine's one contraction entry point. It
does not just multiply the two networks together — it *contracts* them, meaning it multiplies **and
sums over every shared index**:

1. It closes the `dirac` chain into $4\times4$ gamma products and takes the spinor trace (part a's
   `matmul`/`mtrace`, done for you). Legs left open with `dgamma` survive as *free Lorentz indices*
   on the resulting tensor, and the $\mathrm{tr}\,\mathbb{1} = 4$ rides along.
2. It contracts that tensor against the `lorentz` network by index elimination — summing away each
   Lorentz index shared between a `dgamma` leg and a network factor.
3. It returns the surviving scalar as one `MPoly`.

The four arguments: `dirac` = the closed gamma chain; `lorentz` = the `NNet` that ties off the free
legs (empty `{}` when the chain is already scalar); `comp` = the component table; `atomDen` = the
inverse-atom denominators ([step-04](step-04.md)).
```

`Tutorials/step-03-dirac-traces/trace_tokens.cpp`. The frame setup is the same eight symbols:

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
:start-after: "@snip begin: comp"
:end-before: "@snip end: comp"
```

### A chain with no free legs

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
:start-after: "@snip begin: chain-pq"
:end-before: "@snip end: chain-pq"
```

A `DiracNet` is a list of factors **in trace order**, and it is implicitly cyclic — the last
spinor index is tied back to the first, because that is what closing a quark loop means. So this
list of two `dslash` tokens already *is* $\mathrm{tr}(\slashed p\slashed q)$; the `lorentz`
argument is empty because there is nothing left open.

`dslash({{1.0, 0}})` reads "coefficient 1 times momentum id 0". The list is a linear combination,
which is how an internal propagator momentum $q = l - p$ is written: `{{1.0, l}, {-1.0, p}}`. No
new momentum id is introduced for it — the frame only knows the independent momenta.

### A chain with free legs

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
:start-after: "@snip begin: chain-g"
:end-before: "@snip end: chain-g"
```

`dgamma(mu)` is a $\gamma^\mu$ whose Lorentz index is **open**. A chain containing open legs is not
a scalar — it is a tensor with those free indices — so it must be handed a Lorentz network that
ties them off. Here that network is a single metric, which contracts $\mu$ with $\nu$ and gives
$\mathrm{tr}(\gamma^\mu\slashed p\,\gamma_\mu\slashed q)$.

This is the shape of every real diagram: the Dirac chain carries the gamma structure with open
gluon legs, and the Lorentz network carries the gluon propagators and projectors that close them.
[step-05](step-05.md) assembles exactly that.

## Results, part b

```bash
cmake --build build --target trace_tokens && ./build/trace_tokens
```

```text
tr(p/ q/)           = 1.04   (= 4 p.q = 1.04, 4 monomials)
tr(g^mu p/ g_mu q/) = -2.08   (= -2 tr(p/ q/) = -2.08)
ALL TESTS PASSED
```

The first line reproduces part a exactly, as it must — the token path and the matrix path are the
same computation with different bookkeeping.

The second is the 4-dimensional identity $\gamma^\mu\slashed a\,\gamma_\mu = -2\slashed a$, giving
$-2\,\mathrm{tr}(\slashed p\slashed q) = -8\,p\cdot q$. The factor $-2$ is where the Euclidean,
strictly-4-dimensional convention shows: in $d$ dimensions it would be $(2-d)$, and there is no
general-$d$ mode in NumTracer. If your calculation needs dimensional regularisation, this is the
boundary — see [Scope & conventions](../getting_started/scope-and-conventions.md).

Note finally that `numeric_value` returned a *polynomial*, and `eval` put numbers in afterwards.
That separation is the whole basis of code generation: the polynomial is computed once at build
time, and the kernel that ships evaluates it millions of times.

## Possibilities for extensions

1. **Odd numbers of gammas.** Trace a chain of three slashes. You should get exactly zero — every
   trace of an odd number of gammas vanishes. Confirm the polynomial has *no* monomials at all
   (`tr.size() == 0`), not merely monomials that evaluate to zero.

2. **The four-gamma identity.** Verify
   $\mathrm{tr}(\slashed a\slashed b\slashed c\slashed d) = 4[(a\!\cdot\!b)(c\!\cdot\!d) -
   (a\!\cdot\!c)(b\!\cdot\!d) + (a\!\cdot\!d)(b\!\cdot\!c)]$ with four independent momenta
   (`nsym = 16`). Count the monomials before you run it and see whether your guess was right.

3. **$\gamma_5$.** Add `dg5()` to a chain of four slashes. The result is proportional to the
   Levi-Civita tensor and vanishes unless all four momenta are linearly independent — so it is
   *zero* in the 1-angle frames of [step-04](step-04.md) and nonzero with four generic momenta.
   Check both.

4. **Break it deliberately — forget to close the legs.** Build `chainG` but pass `{}` as the
   Lorentz network. Read the failure carefully: an open Lorentz index with nothing to contract it
   against is a real error and the engine should tell you so. Compare with what happens if you pass
   a metric on the *wrong* labels — that one is not an error at all, just a different (and silently
   wrong) contraction.

5. **Measure the scaling.** Time chains of 4, 8, 12, 16 slashes. The matrix-product cost grows
   linearly in the chain length; the Wick pairing count grows as $(2n-1)!!$. Plot both on a log
   axis and you have the argument for the whole design.

## The plain programs

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.cpp
:language: cpp
```

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
```

Next, [step-04](step-04.md): the Lorentz network that ties those open legs off, and what happens to
the $1/l^2$ it drags along.
