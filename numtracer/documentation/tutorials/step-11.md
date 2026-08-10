# step-11: Projectors and tensor bases

*Builds on: [step-08](step-08.md), [step-10](step-10.md) · Built on by: [step-12](step-12.md),
[step-14](step-14.md) · Tags: `tensorbases`, `projector` · **Tier C** (the full fRG toolchain)*

````{admonition} What this step needs
:class: important
FunKit, TensorBases, DiFfRG and FORM. Not wired into `ctest`; run it as

```bash
NT_TUTORIAL_OUT=/tmp/out wolframscript -script Tutorials/step-11-projectors/orthonormality.wl
```
````

## Introduction

[step-10](step-10.md) projected the quark two-point flow with this line:

```mathematica
FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRGqbq["1-Loop"]["ExternalIndices"]]] ** …
```

and three conventions were asserted in passing: a *minus sign*, a *reversed index order*, and the
choice of `TBGetProjector` over its neighbours. **NumTracer cannot check any of them.** It contracts
whatever network you hand it. If the projector is transposed, the kernel is still a valid kernel; it
just computes a different number, and you find out months later when a coupling runs the wrong way.

This step is about the check that *does* catch it — and it needs no tracer, no flow, and no oracle.

### Why a projector exists at all

A vertex is a tensor. The flow of a three-gluon vertex is a rank-3 Lorentz tensor with colour
indices; the flow of a quark two-point function is a matrix in Dirac ⊗ colour ⊗ flavour space. But a
kernel returns a *number*, and a numerical fRG solves for a handful of *scalar* dressings.

The bridge is a **tensor basis**. You choose a set of structures $T_i$ spanning the tensors the
symmetry allows, and write

$$
\Gamma^{(n)} \;=\; \sum_i \lambda_i(p)\; T_i ,
$$

so the physics lives in the scalar coefficients $\lambda_i$. To extract one you need a **dual
basis** $P_i$ — the projectors — defined by exactly one property:

$$
\langle P_i,\, T_j\rangle \;=\; \delta_{ij}.
$$

Constructing $P_i$ means inverting the Gram matrix $\langle T_i, T_j\rangle$, which TensorBases does
(with FORM) and caches under `TBCache/`.

That defining property is the whole check. If it holds, your projectors are the duals of your
structures and every extracted $\lambda_i$ means what you think. If it fails, nothing downstream is
salvageable.

### The accessors, and which one to use

| Accessor | Returns |
|---|---|
| `TBGetBasisElement[b, i, idx]` | the structure $T_i$ itself |
| `TBGetVertex[b, i, idx]` | the *vertex* built from element $i$ |
| `TBGetProjector[b, i, idx]` | the dual $P_i$ |
| `TBGetMetric` / `TBGetInverseMetric` | the Gram matrix and its inverse |
| `TBGetInnerProduct[b, x, y]` | the basis's inner product |

```{admonition} A projector is dual to the VERTEX
:class: warning
Pair `TBGetProjector` with **`TBGetVertex`**, and contract them directly with indices you supply by
hand. Do **not** route the check through `TBGetInnerProduct`: for a vertex basis (`AqbqDirect` is
one) `TBMakeInnerProduct` selects the canonical product used to project *onto vertices*, so pairing
it with `TBGetBasisElement` mixes two conventions and gives a plausible-looking wrong answer.
```

## The commented program

`Tutorials/step-11-projectors/orthonormality.wl`, on the two-element `"qbq"` basis — the same one
[step-10](step-10.md) used.

### First: guard against a vacuous pass

```{literalinclude} ../../../Tutorials/step-11-projectors/orthonormality.wl
:language: mathematica
:start-after: "@snip begin: guard"
:end-before: "@snip end: guard"
```

This is the most important part of the script. If a basis never resolved — no cached Gram matrix and
FORM did not build one — `TBGetProjector` stays **unevaluated**. It is then just a symbolic head,
and a symbolic head flows through `FromFunKit`, through `NumTrace`, into a kernel, without ever
being a projector. Every check below would "pass" while testing nothing.

A green orthonormality gate on an unresolved basis is worse than no gate at all.

Note also the index bundles:

```mathematica
idxV = {{p1, d1, A1, F1}, {-p1, d2, A2, F2}};
idxP = Reverse[idxV];
```

Each entry is one leg's **complete bundle**: momentum, Dirac, colour, flavour. Momentum
conservation forces the second leg to carry $-p_1$. And the momentum **travels with its leg** — you
move whole bundles, never just the index triples.

### The check

```{literalinclude} ../../../Tutorials/step-11-projectors/orthonormality.wl
:language: mathematica
:start-after: "@snip begin: pairs"
:end-before: "@snip end: pairs"
```

`FormTrace` closes the Dirac, colour and flavour indices. The result must be the identity matrix —
**exactly**, symbolically, at generic momenta. Not "to $10^{-12}$": these are rational functions,
and a correct dual basis gives literal `1` and `0`.

### What the wrong convention looks like

```{literalinclude} ../../../Tutorials/step-11-projectors/orthonormality.wl
:language: mathematica
:start-after: "@snip begin: wrong"
:end-before: "@snip end: wrong"
```

Rather than assert the rule, the script *demonstrates* it: it repeats the check with the projector
given the **same** leg order as the vertex — the natural-looking thing to write — and prints what
comes back.

### Scaling up: the weighted-trace trick

```{literalinclude} ../../../Tutorials/step-11-projectors/orthonormality.wl
:language: mathematica
:start-after: "@snip begin: weighted"
:end-before: "@snip end: weighted"
```

A basis with $n$ structures has $n^2$ pairs. For the full 12-element quark–gluon vertex basis that
is 144 — and if you want the check to run through the *generated kernel* rather than in Mathematica,
144 committed headers. Instead, contract the whole matrix against a fixed weight matrix in one net:

$$
W \;=\; \sum_{ij} w_{ij}\,\langle P_i, V_j\rangle \quad\text{must equal}\quad \sum_i w_{ii}.
$$

The weights are **distinct primes**, so no plausible error cancels: any single wrong entry shifts
$W$ by its own weight, and the diagonal and off-diagonal weights are disjoint. One kernel covers all
$n^2$ pairs. The per-pair matrix is then a *diagnostic* — what you run when $W$ is wrong and you
need to know which entry moved.

This is exactly what `numtracer/tests/gen/gen_ortho_numeric.wls` does for the real bases, and it is
worth reading once you have this small version straight.

## Results

```text
step-11: basis "qbq" resolved.
step-11: <P_i, V_j> = {{1, 0}, {0, 1}}
step-11: orthonormality holds exactly.
step-11: with the projector's legs NOT reversed, <P_i, V_j> =
         {{-1/4*gamma[p1, spi1, spi2]^2/sp[p1, p1], 0}, {0, 1}}
step-11: weighted trace = 9   (expected 9)
step-11: all projector checks passed.
```

Read the third line carefully, because recognising this shape is worth more than memorising the
rule. With the legs not reversed, the $(1,1)$ entry came back as

$$
-\frac{1}{4}\,\frac{\gamma[p_1, \mathrm{spi}_1, \mathrm{spi}_2]^2}{p_1^2}
$$

— an **unclosed spinor line**. The two gamma matrices never paired into a trace, so `FormTrace` left
a squared open-index object standing. That is the signature of a transposed Dirac structure.

Note what *did* still work: the $(2,2)$ entry is 1 in both conventions, because the mass structure is
proportional to the identity in Dirac space and does not care about transposition. **A partially
correct answer is the normal appearance of this bug** — which is why the check must cover every
pair, not just the one you happen to be extracting.

## Two conventions worth knowing before they bite

### The momentum-routing sign

`QCD_Nf2/baseline/QCD.wl` carries a 20-line comment that is the best single illustration of how far
a projector convention propagates:

```{admonition} From QCD_Nf2/baseline/QCD.wl — a basis-element parity that reversed the quark loop
:class: warning
> TensorBases 1.3.0 (commit 9fbdc0a, "Make TBMakePropagator use all-incoming momenta") changed how
> the InverseProp dressings below are inverted … Element 1 of the `"qbq"` basis is
> `I*pdash[p1,d1,d2]` — **MOMENTUM-ODD** — so its returned propagator dressing flips sign under the
> new routing, while the momentum-even mass element 2 and every bosonic basis are untouched.
> Regenerating against 1.3.0 without compensating therefore flips every trace term with an odd
> number of quark pslash insertions … which reverses the quark loop everywhere: `m2A` then runs UP
> instead of toward criticality for any initial value, `etaQ` flips, and `etaPiL` flips so that
> $Z_\pi$ COLLAPSES instead of growing — the flow dies at $t \approx 0.85$.

An upstream change in *momentum routing*, interacting with the *parity of one basis element*,
reversed the sign of a whole class of trace terms and killed the flow. The compensating minus lives
in `dressingRules` and nowhere else. Note that the failure was in the physics, not in any test.
```

### Asserting a hand-built projector

Not every projection comes from a basis. The Yukawa flow in `QCD_Nf2/with_mesons` builds its
projector by hand — and then immediately checks it:

```mathematica
hPiStruct    = TFlav[f1, F2, F3]*((-I)*gamma5[d2, d3])*deltaFundCol[A2, A3];
hPiStructBar = TFlav[f1, F3, F2]*((+I)*gamma5[d3, d2])*deltaFundCol[A3, A2];
NormhPi      = FormTrace[hPiStruct*hPiStructBar];
ProjectorhPi = hPiStructBar/NormhPi;

Module[{sanity},
  sanity = Simplify[FormTrace[ProjectorhPi *
     ((GammaN[{ΠL, qb, q}, {…}] /. diagRules) /. hPhiL[__] :> 1)]];
  Print["ProjectorhPi . Gamma^(3)/h = ", sanity, "  -- must be 1"];
  If[sanity =!= 1, Print["!! Yukawa projector normalisation/sign is WRONG -- aborting"]; Abort[]]];
```

Two things to copy. The dual **transposes every index pair and conjugates**, so that the Dirac pair
closes into a trace — `gamma5[d2,d3]*gamma5[d2,d3]` does *not* close, and FORM leaves it as an
unclosed square (the same signature as above). And the check `Abort[]`s. A printed warning nobody
reads is not a check.

The same file's history shows what the alternative costs:

```{admonition} From QCD_Nf2/with_mesons_taylor/QCD.wl
:class: warning
> The LEADING MINUS is not cosmetic. … Without the minus here,
> `ProjectorhPhiL . Gamma^(2) = -hPhiL` … so the whole quark self-energy — i.e. the gluon-exchange
> generation of the quark mass, THE driver of chiral symmetry breaking — entered `dt h` with the
> wrong sign and DRAINED the meson channel instead of feeding it.
```

### Degenerate bases and the symmetric point

`gen_ortho_numeric.wls` uses `gen3Frame`, not `sp3Frame`, and says why:

> `sp3Frame` pins a measure-zero slice on which the full `AqbqDirect` basis is degenerate
> (`Det[Gram] = 0` there, so every dual projector is 0/0); checking orthonormality only at that
> point would test the tracer where the algebra breaks down. Verified: `Det[Gram] = 0` at
> $(1,1,-1/2)$ and $1.87\times10^{22}$ at $(1, 7/4, -3/10)$.

This is the concrete reason [step-07](step-07.md)'s `gen3Frame` exists. If your basis has a singular
Gram at the symmetric point, the symmetric point is not a valid place to work — and an
orthonormality check performed there will not tell you so.

## Possibilities for extensions

1. **Break the guard.** Point `SetDirectory` somewhere without a `TBCache/` and disable FORM. Watch
   `TBGetProjector` come back unevaluated, then comment out the guard and confirm the rest of the
   script "passes". This is the failure the guard exists for.

2. **Transpose only the indices.** Swap the `(d, A, F)` triples between the two bundles while
   pinning the momenta — i.e. violate "the momentum travels with its leg". Compare the resulting
   matrix with the reversed-legs case. On this basis they differ; on a symmetric bosonic basis they
   may not, which is why the script asserts its own demonstration is non-vacuous.

3. **Do it on a real vertex basis.** Repeat with `"AqbqDirect147"` (3 structures, 9 pairs) using
   `gen3Frame` kinematics. Then check `Det` of the Gram at the symmetric point versus a generic
   point and reproduce the numbers quoted above.

4. **Push it through the kernel.** The check above runs in Mathematica. Wire the weighted trace
   through `NumTrace` + `MakeNTKernel` and grade the emitted kernel at random momenta, as
   `gen_ortho_numeric.wls` does. That upgrades it from "the basis is consistent" to "the basis is
   consistent *and the tracer contracts it correctly*" — a genuinely stronger statement.

5. **Reproduce the sign disaster.** In a *copy* of `qcd_setup.m`, restore the minus on the
   `InverseProp {qb,q} 1` rule. Regenerate a `Zq` kernel and compare against the committed one. Then
   consider that no test in the repository would have caught this — only the physics did.

## The plain program

```{literalinclude} ../../../Tutorials/step-11-projectors/orthonormality.wl
:language: mathematica
```

Next, [step-12](step-12.md): a complete propagator flow — the ghost and gluon two-point functions —
where the projector, the symmetry reduction and the dressing rules all have to be right at once.
