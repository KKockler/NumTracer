# step-12: A propagator flow, end to end

*Builds on: [step-10](step-10.md), [step-11](step-11.md) · Built on by: [step-13](step-13.md) ·
Tags: `frg`, `flow` · **Tier C** (the full fRG toolchain)*

```{admonition} What this step needs
:class: important
FunKit, TensorBases, DiFfRG and FORM. Rather than ship a copy, this step reads the **production**
code: `QCD_Nf2/baseline/QCD.wl` — a complete two-flavour QCD fRG project in 588 lines. Every listing
below is quoted from it verbatim.
```

## Introduction

[step-10](step-10.md) derived one flow. This step reads two more — the **ghost** and **gluon**
two-point functions — from real production code, because they are where three things that have been
introduced separately finally have to be right at the same time: the projector and its sign
([step-11](step-11.md)), the leg symmetries, and the dressing rules.

They are also, deliberately, almost the same code. The value of reading them side by side is in the
differences.

### The shape every flow has

Every one of the 7 flows in `baseline`, and all 31 in `with_mesons`, is the same five lines:

```mathematica
Module[{sym, fRG, tr, ntk},
  sym = FMakeSymmetryList[...];                                  (* 1. leg symmetries, optional *)
  fRG = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]},
          "Symmetries" -> sym] // FTruncate // FPlot // FRoute;  (* 2. differentiate *)
  tr  = FTerm[prefactor, TBGetProjector[...],
              (fRG["1-Loop"]["Expression"] /. diagRules)];       (* 3. project *)
  ntk = NumTrace[FromFunKit[PropParam[dressingRules[tr]]],
                 "Frame" -> frame1, "Args" -> args1];            (* 4. trace *)
  MakeNTKernelDiFfRG[ntk, "Name" -> "ZA", ...];                  (* 5. emit into DiFfRG *)
];
```

Learn this shape and you can read any flow in the tree. Steps 10 and 11 covered stages 2–4; stage 5
is [step-15](step-15.md).

## The gluon two-point flow

```mathematica
Module[{sym, fRG, tr, ntk},
  (* the transverse gluon 2-point projector is symmetric in the two gluon legs *)
  sym = FMakeSymmetryList[FSymmetry[Symmetric, {i1, i2}]];
  fRG = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}, "Symmetries" -> sym]
          // FTruncate // FPlot // FRoute // FPrint;
  tr  = FTerm[1/sp[p1, p1],
              TBGetProjector["AA", 1, {i1, i2} /. fRG["1-Loop"]["ExternalIndices"]],
              (fRG["1-Loop"]["Expression"] /. diagRules)] //. p2 -> -p1;
  ntk = NumTrace[FromFunKit[PropParam[dressingRules[tr]]], "Frame" -> frame1, "Args" -> args1];
  Print["ZA nets = ", Length[ntk[[1]]["Diagrams"]]];
  MakeNTKernelDiFfRG[ntk,
      "Name" -> "ZA", "Namespace" -> "za_qcd",
      "Integrator" -> "Integrator_p2_1ang",
      "Parameters" -> kernelParameterList,
      "IntegrationVariables" -> {"l1", "cos1"}];
];
```

Four things to notice.

**`//. p2 -> -p1`.** Momentum conservation, applied *after* the projection. A two-point function has
one independent external momentum.

**`FPlot` and `FPrint`** in the pipeline. They draw the diagrams and print the expression. Neither
changes anything — they are `Identity` with a side effect — and they are in the production code
because *looking at the diagrams* is how you find out that your truncation is not what you meant.
Use them.

**The `1/sp[p1, p1]` prefactor with a `+` sign.** This is the normalisation of the gluon dressing:
`dressing[GammaN, {A,A}, 1, …] :> ZA[…] sp[p2,p2]` carries a factor of $p^2$, so the projection
divides it out. Hold that thought.

**`"Integrator" -> "Integrator_p2_1ang"`** matches `frame1 = propFrame[…]` — one loop magnitude,
one angle. The frame and the integrator must agree about how many integration variables there are;
nothing checks it ([step-15](step-15.md)).

### When is the symmetry reduction legitimate?

The `"Symmetries"` option is not free. The production file carries a careful argument for why it is
allowed, and it deserves to be read in full:

```{admonition} From QCD_Nf2/baseline/QCD.wl
:class: important
> **SYMMETRIES.** The reduction happens inside `FTruncate`, on the bare diagram sum with symbolic
> external momenta — BEFORE the projector is attached and before `dressingRules`. It replaces
> $D + \sigma(D)$ by $2D$, so the reduced sum **differs from the exact one as an expression**; it is
> usable iff everything applied afterwards annihilates the difference, i.e. iff
> $P\cdot\sigma(X) = P\cdot X$ for all $X$. That is a property of the projector as a tensor function
> of indices AND momenta at GENERIC momenta — the kinematic point (`p2 -> -p1`, the SP3/SP4 frames)
> is then irrelevant, since it is substituted into both sides of an identity. The second condition
> is that `dressingRules`, applied after the reduction, be covariant too: it is, because every
> vertex parametrization here is a symmetric momentum average (`SP3Patt`, `SP4Patt`) and every
> propagator dressing depends only on its own leg. **An asymmetrically parametrized vertex added
> later would silently invalidate all of this.**
```

Three things worth extracting:

1. The reduction is an *expression-level* rewrite that is only valid *given what you do next*. It is
   not a symmetry of the object; it is a symmetry of the object-as-projected.
2. The condition is at **generic momenta**. Checking it at your kinematic point proves nothing,
   because the point is substituted into both sides of an identity that may not hold.
3. The last sentence is the real warning. This is a correctness condition that a *future* change —
   adding a vertex with an asymmetric parametrisation — silently breaks, in a file that will not be
   re-read. Note that no test catches it.

There are two spellings of the option, and they are not equivalent:

| Spelling | Group generated |
|---|---|
| `FMakeSymmetryList[FSymmetry[Symmetric, {i1, i2}]]` | exactly the stated permutation |
| `FMakeSymmetryList[{A[i1], A[i2], A[i3]}]` (field list) | the **full** permutation group of those legs |

`ZA3` uses the field-list form, because a three-gluon vertex basis is Bose-symmetric in all three
legs at generic momenta. [step-19](step-19.md) shows a case where the field-list form is *wrong* and
costs a factor of five in respected diagrams.

## The ghost two-point flow

Almost the same code. The differences are the lesson:

```mathematica
Module[{fRG, tr, ntk},
  (* no symmetry: cb and c are distinct fields, there is no leg permutation to exploit *)
  fRG = FTakeDerivatives[WetterichEquation, {cb[i1], c[i2]}] // FTruncate // FPlot // FRoute // FPrint;
  (* NOTE the MINUS on 1/sp: the ghost dressing is defined with one, i.e.
     dressing[GammaN, {cb, c}, 1, ...] :> -Zc[|p|]*sp[p2, p2], so the projection has to
     carry the compensating -1/sp[p1, p1]. (ZA has +ZA[|p|]*sp and hence +1/sp.)
     Dropping it flips the sign of the whole ghost flow. *)
  tr  = FTerm[-1/sp[p1, p1],
              TBGetProjector["cbc", 1, {i1, i2} /. fRG["1-Loop"]["ExternalIndices"]],
              (fRG["1-Loop"]["Expression"] /. diagRules)] //. p2 -> -p1;
  ntk = NumTrace[FromFunKit[PropParam[dressingRules[tr]]], "Frame" -> frame1, "Args" -> args1];
  MakeNTKernelDiFfRG[ntk, "Name" -> "Zc", "Namespace" -> "zc_qcd", …];
];
```

* **No symmetry.** $\bar c$ and $c$ are distinct fields; there is no leg permutation.
* **The prefactor is $-1/p^2$, not $+1/p^2$.** Not because the physics differs, but because
  `dressingRules` defines the ghost dressing with a minus:
  `dressing[GammaN, {cb,c}, 1, …] :> -Zc[…] sp[p2,p2]`. The projection must compensate.

That minus sign is a pure bookkeeping convention with no physical content, and getting it wrong
flips the entire ghost flow. It is invisible to every check in the toolchain — the kernel compiles,
the orthonormality gate of [step-11](step-11.md) passes (it tests the *basis*, not your prefactor),
and the flow simply runs the wrong way.

**The general rule:** your projection prefactor must undo exactly the normalisation your
`dressingRules` put on that dressing. Write them next to each other and check the pair, every time.

## The quark two-point flow, and why it has no prefactor

`Zq` is where the pattern breaks, and the production comment explains the consequence precisely:

```{admonition} From QCD_Nf2/baseline/QCD.wl
:class: warning
> **NOTE: NO `1/sp` prefactor here, unlike ZA/Zc.** The quark dressing is
> `dressing[GammaN, {qb, q}, 1, ...] :> -zq[|p|]` — **no factor of `sp[p2, p2]`** — so an extra
> `1/sp[p1, p1]` would over-normalise the flow by $1/p^2$. At `p_grid_min = 1e-3` that is a **1e6
> enhancement** at the pinning point, which turns the `etaQ` Newton iteration in `model.hh` into an
> expanding map (`etaQ -> 1e250 -> NaN` on the first timestep).
```

This is the same rule as the ghost's, and the failure mode is spectacular rather than subtle — a
NaN on the first timestep instead of a slowly wrong answer. Both come from the same question:
*what normalisation did `dressingRules` give this dressing?*

The minus that `ZA`/`Zc` carry on the prefactor sits on the **projector** instead here, together
with the reversed index order `{i2, i1}` that [step-11](step-11.md) demonstrated:

```mathematica
tr = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRG["1-Loop"]["ExternalIndices"]],
           (fRG["1-Loop"]["Expression"] /. diagRules)];
```

The positional second argument to `MakeNTKernelDiFfRG` — the counterterm `etaQ*zq[p]` — is
[step-13](step-13.md).

## Reading the emitted kernel

`no_mesons/flows/ZA/kernels.hh` has `nenv = 4` and **eight** trace functions, and the kernel's
return is a single nested `fma` chain over them. Stripped of the CSE variables:

```cpp
return fma(-0.0833…, … tr3(fenv) * _interp13,
       fma(-0.1666…, … tr4(fenv) * powr<2>(_interp15),
       fma(-0.1666…, … tr5(fenv) * powr<2>(_interp18),
       fma(-0.0833…, … tr6(fenv) * _interp15 * _interp19,
       fma(-0.0416…, … tr7(fenv) * powr<2>(_interp19),
       fma(-0.0416…, … tr0(fenv) * powr<2>(_interp6),
       fma(-0.0208…, … tr1(fenv) * _interp7,
       fma(-0.0833…, … tr2(fenv) * powr<2>(_interp8)
                       * (_interp2*nthk3 + _interp1*(nthk4 + 50.*(-nthk3 + nthk5))),
       0.))))))));
```

**Each `trN(fenv)` is one diagram** (or one group of diagrams sharing a coefficient). You are looking
at the gluon self-energy: a gluon tadpole from the four-gluon vertex, a ghost loop, a gluon loop, a
quark loop, and — in `no_mesons` — the two-quark-two-gluon contributions. The
$\sum_{\text{diagrams}}\text{coeff}\times\text{trace}$ shape from [step-06](step-06.md) is still
exactly what you see, now with eight terms.

**The `_interpN` are dressing evaluations**, hoisted and computed once. **The `nthkN`** are the
finite-difference stencil points the regulator's $\partial_t$ needs.

**`nenv = 4` for eight traces.** All eight share one four-slot environment: the frame symbols are
filled once per call and every trace reads them. That sharing is why adding a diagram to a flow
costs much less than the first one did.

## Possibilities for extensions

1. **Diff `Zc` against `ZA`.** Put the two `Module[…]` blocks side by side. There are exactly three
   differences: the symmetry list, the sign of the prefactor, and the basis name. Every other
   character is identical. That is what "the shape every flow has" means in practice.

2. **Flip the ghost sign.** In a copy of the tree, change `-1/sp[p1,p1]` to `+1/sp[p1,p1]`,
   regenerate, and run the flow. It will not crash. Watch $Z_c$ run the wrong way and think about
   how you would have noticed if you had not been told.

3. **Check a normalisation pair.** For every dressing in `dressingRules`, write down the power of
   $p^2$ it carries and the prefactor the corresponding flow uses. They should be reciprocal. Any
   pair that is not is either a bug or a convention worth a comment.

4. **Test the symmetry condition.** For the `AA` projector, verify $P\cdot\sigma(X) = P\cdot X$
   symbolically at generic momenta for a general tensor $X$. Then do it for a vertex with an
   asymmetric momentum parametrisation and watch it fail — the failure the production comment warns
   a future maintainer about.

5. **Count the diagrams.** `Print["ZA nets = ", Length[ntk[[1]]["Diagrams"]]]` is in the production
   code for a reason. Run `baseline` and `no_mesons` and compare: the extra quark-gluon tensors and
   the two-quark-two-gluon vertex are visible as a jump in the count long before anything is
   emitted.

Next, [step-13](step-13.md): the counterterm in `Zq`'s positional slot, and the Newton iteration in
`model.hh` that cannot converge without it.
