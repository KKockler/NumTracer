# step-14: Three- and four-point vertices

*Builds on: [step-11](step-11.md), [step-13](step-13.md) · Built on by: [step-15](step-15.md),
[step-19](step-19.md) · Tags: `vertex`, `frame` · **Tier C** (the full fRG toolchain)*

```{admonition} What this step needs
:class: important
Read from `QCD_Nf2/baseline/QCD.wl` and `QCD_Nf2/no_mesons/`. All listings are quoted verbatim.
```

## Introduction

Every flow so far has been a two-point function: one external momentum, one loop angle, a handful
of diagrams. Vertices are where the calculation becomes expensive, and where three mechanisms that
have been mentioned in passing become essential.

The scale, measured on `no_mesons/flows/*/kernels.hh`:

| Flow | External legs | Loop angles | Emitted size | Trace functions |
|---|---|---|---|---|
| `Zc` | 2 | 1 | 1.7 kB | 2 |
| `ZA` | 2 | 1 | 4.1 kB | 8 |
| `ZAcbc` | 3 | 2 | 11 kB | 6 |
| `ZA3` | 3 | 2 | 39 kB | 3 |
| `ZAqbq1` | 3 | 2 | 124 kB | 55 |
| `ZAqbq4` | 3 | 2 | 194 kB | 102 |
| `ZAAqbq1` | 4 | 3 | 59 kB | 14 |
| `ZA4` | 4 | 3 | 467 kB | 75 |

Two orders of magnitude between the first and last row. In `with_mesons` the same `ZA4` reaches
5 MB and `lambda1L3D` reaches **81 MB** with 989 trace functions. Everything in
[step-20](step-20.md) exists because of this table.

## The simplest vertex: three distinct fields

`ZAcbc` — a gluon and a ghost pair. Structurally the easiest possible vertex, because three
distinct fields admit no leg permutation at all:

```mathematica
Module[{fRG, tr, ntk},
  (* no symmetry: A, cb, c are three distinct fields *)
  fRG = FTakeDerivatives[WetterichEquation, {A[i1], cb[i2], c[i3]}] // FTruncate // FPlot // FRoute;
  tr  = FTerm[TBGetProjector["Acbc", 1, {i1, i2, i3} /. fRG["1-Loop"]["ExternalIndices"]],
              (fRG["1-Loop"]["Expression"] /. diagRules)];
  ntk = NumTrace[FromFunKit[SPParam[dressingRules[tr]]], "Frame" -> frame2, "Args" -> args2];
  MakeNTKernelDiFfRG[ntk,
      "Name" -> "ZAcbc", "Namespace" -> "zacbc_qcd",
      "Integrator" -> "Integrator_p2_4D_2ang",
      "Parameters" -> kernelParameterList,
      "IntegrationVariables" -> {"l1", "cos1", "cos2"},
      "AngleDefs" -> spAngles3];
];
```

Compare against `ZA` from [step-12](step-12.md). Four things changed, and they change together:

| | two-point | three-point |
|---|---|---|
| frame | `frame1` = `propFrame` | `frame2` = `sp3Frame` |
| args | `{l1, cos1, p, k}` | `{l1, cos1, cos2, p, k}` |
| scalar reduction | `PropParam` | `SPParam` |
| integrator | `Integrator_p2_1ang` | `Integrator_p2_4D_2ang` |
| angle map | — | `"AngleDefs" -> spAngles3` |
| prefactor | `±1/sp[p1,p1]` | none |

**These five must agree with each other, and nothing checks that they do.** A `frame2` with an
`Integrator_p2_1ang`, or `SPParam` with `spAngles4`, produces a kernel that compiles and integrates
the wrong thing.

### `SPParam` and `AngleDefs`: the two halves of the kinematics

This is the mechanism [step-07](step-07.md) introduced and here it does real work.

The **scalar** half is reduced by `SPParam`, which — unlike `PropParam` — leaves the loop–external
cosines *symbolic*:

```mathematica
SPParam[expr_] := UseLorentzLinearity[expr] //. {
  lf1 -> l1, sp[p, p] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> p l1 cos[l1, p1], sp[l1, p2] -> p l1 cos[l1, p2], …
  cos[l1, p1] :> cosl1p1, cos[l1, p2] :> cosl1p2, cos[l1, p3] :> cosl1p3};
```

The **tensor** half gets its kinematics from `frame2` instead. `"AngleDefs"` is the bridge:

```mathematica
spAngles3 = {
  cosl1p1 -> cos1,
  cosl1p2 -> (-cos1 + Sqrt[3 - 3 cos1^2] cos2)/2,
  cosl1p3 -> (-cos1 - Sqrt[3 - 3 cos1^2] cos2)/2};
```

Those are exactly $\hat l\cdot\hat q_i$ computed from the `sp3Frame` component table printed in
[step-07](step-07.md) — derive them and check. Because they are *derived from the frame*, an
`AngleDefs` written for one frame and used with another silently puts the two halves of the kernel
on different kinematics.

## Bose symmetry: the field-list form

`ZA3` differs from `ZAcbc` in exactly one line:

```mathematica
(* AAAClassTrans is Bose-symmetric in all three gluon legs (indices and momenta together,
   at generic momenta), so the FULL permutation group of the three legs may be used *)
sym = FMakeSymmetryList[{A[i1], A[i2], A[i3]}];
```

The **field-list form** generates the full permutation group of those legs — six permutations here,
versus the single transposition `FSymmetry[Symmetric, {i1, i2}]` gave `ZA` in
[step-12](step-12.md). The saving is real: `ZA3` emits 3 trace functions where `ZAcbc` needs 6,
despite being the larger object.

The parenthetical is the licence. The permutation must be a symmetry **of the indices and momenta
together, at generic momenta** — the condition [step-12](step-12.md) quoted in full. Do not reach
for the field-list form because the legs carry the same field; reach for it because you have checked
that the projector is invariant. [step-19](step-19.md) shows a four-fermi case where the same-field
reasoning gives the wrong group and loses 82 % of the diagrams.

## One derivative, several projectors

The quark–gluon vertex has three tensor structures in play (1, 4, 7 of the `AqbqDirect` basis). All
three share the *same* functional derivative, so `baseline/QCD.wl` memoises it:

```mathematica
(* no symmetry: A, qb, q are three distinct fields *)
zaqbqFRG := zaqbqFRG = FTakeDerivatives[WetterichEquation, {A[i1], qb[i2], q[i3]}]
              // FTruncate // FPlot // FRoute;
```

`:=` with the body assigning to itself is Mathematica's memoisation idiom: the derivative is taken
once, on first use, and cached. `no_mesons` then emits three kernels off it:

```mathematica
Module[{fRG, tr, ntk},
  fRG = zaqbqFRG;
  tr  = FTerm[TBGetProjector["AqbqDirect", 4, {i1, i3, i2} /. fRG["1-Loop"]["ExternalIndices"]],
              (fRG["1-Loop"]["Expression"] /. diagRules)];
  ntk = NumTrace[FromFunKit[SPParam[dressingRules[tr]]], "Frame" -> frame2, "Args" -> args2];
  MakeNTKernelDiFfRG[ntk, etaQ*ZAqbq4[p], "Name" -> "ZAqbq4", …];
];
```

**Note the projector index order `{i1, i3, i2}`** — the quark legs are swapped relative to the
derivative order `{A[i1], qb[i2], q[i3]}`, the three-point version of the `{i2, i1}` reversal
[step-11](step-11.md) demonstrated.

### Restricting the basis

Which structures exist at all is set once, in the setup:

```mathematica
TBRestrictBasis["AqbqDirect", "AqbqDirect1",   {1}];        (* baseline: classical tensor only *)
TBRestrictBasis["AqbqDirect", "AqbqDirect147", {1, 4, 7}];  (* no_mesons: three structures *)
```

`TBRestrictBasis` builds a *new named basis* from a subset of an existing one — and crucially it
recomputes the Gram matrix and the duals for that subset. The projectors of `AqbqDirect1` are **not**
the projectors of `AqbqDirect` with rows deleted; a restricted basis has its own duals. This is why
the truncation ladder in `QCD_Nf2` is a one-line change per tier.

Note that `ZAqbq1`/`ZAqbq4` also carry a counterterm (`etaQ*ZAqbq1[p]`) — a three-point function
with one quark bilinear, exactly as [step-13](step-13.md) predicts.

## Restricting a trace

A basis restriction removes structures from the *external* projection. Restricting what appears
*inside* the loop is a different operation, done with plain rewrite rules:

```mathematica
dropAqbq47[e_] := e //. {ZAqbq4[__] :> 0, ZAqbq7[__] :> 0};

(* Same restriction, but on the RAW dressing heads, i.e. BEFORE dressingRules maps them onto
   ZAqbq4/ZAqbq7. Used for the four-quark box/triangle, whose only consumer is the hadronised
   Yukawa — the hPhiL loop itself keeps the full quark-gluon vertex. *)
dropRawAqbq47[e_] := e //. {dressing[GammaN, {A, qb, q}, 2, _] :> 0,
                            dressing[GammaN, {A, qb, q}, 3, _] :> 0};

dropFourQuark[e_] := e //. {dressing[GammaN, {qb, qb, q, q}, _, _] :> 0};
```

**There are two places to cut**, and the distinction matters: before `dressingRules` you match the
abstract `dressing[GammaN, {A,qb,q}, 2, _]` heads; after it, the named `ZAqbq4[__]`. They are the
same physics and different code, and mixing them up produces a rule that never fires — silently, as
[step-10](step-10.md) warned.

These are applied inline in the pipeline:

```mathematica
ntk = NumTrace[FromFunKit[SPParam[dressingRules[tr] // dropAqbq47 // dropFourQuark]],
               "Frame" -> frame3, "Args" -> args3];
```

The cost is measured. `with_mesons` parameterises the choice by environment variable and records
what each setting costs:

```mathematica
$fqTensors = ToExpression /@ envList["NT_FQ_TENSORS", {"1", "4", "7"}];
```

> 133 / 439 / 1349 diagrams and 13 / 38 / 87 s at ~0.2 GiB for {1} / {1,7} / {1,4,7}.

A tenfold diagram count for one extra tensor structure. This is the lever [step-20](step-20.md)
reaches for first when a flow will not generate.

## Four-point: partial symmetry

`ZAAqbq` has two gluons and a quark pair. The gluons may be exchanged; the quark legs may not:

```mathematica
(* transAAqbqMinimal is symmetric in the two gluon legs; the quark legs are distinct fields
   and admit no permutation, so the field-list form yields exactly {identity, i1<->i2}. *)
symsAAqbq := symsAAqbq = FMakeSymmetryList[{A[i1], A[i2], qb[i3], q[i4]}];
```

Here the field-list form is safe *because* it yields only what you want — the grouping by field type
happens to coincide with the actual symmetry. The comment says so explicitly rather than relying on
the reader to work it out.

Four external legs means `frame3` = `sp4Frame`, three loop angles, `"AngleDefs" -> spAngles4` and
`Integrator_p2_4D_3ang`. The pattern is unchanged; only the numbers grow.

## An aside: a flow with no derivative at all

Not every flow comes from `FTakeDerivatives`. The effective potential is the Wetterich equation
itself, sliced:

```mathematica
emitV[entries_, name_, ns_] := Module[{fRG, tr, ntk},
  fRG = (WetterichEquation // FTruncate)[[entries]] // FPlot // FRoute // FPrint;
  tr  = FTerm[1, (fRG["1-Loop"]["Expression"] //. diagRules)];
  …
  ntk = NumTrace[FromFunKit[PropParam[dressingRules[tr]] //. p1 -> 0], "Frame" -> frame0, …];
  MakeNTKernelDiFfRG[ntk, "Name" -> name, "Namespace" -> ns,
      "Integrator" -> "Integrator_p2",        (* zero angles *)
      "Coordinates" -> {}, "CoordinateArguments" -> {},
      "Device" -> "TBB", "AD" -> True]];

emitV[{3, 4, 5}, "V",    "vfe_qcd"];
emitV[{3, 4},    "Vadv", "vadv_qcd"];
emitV[{5},       "Vsig", "vsig_qcd"];
```

> Supertrace entries are ordered as the `Propagator` list of `truncation`, i.e.
> {A, ghost, quark, pion, sigma}, so {3,4} = {quark, pion} is EXACTLY the advection flux,
> {5} = {sigma} is EXACTLY the diffusion flux, and V == Vadv + Vsig identically.

Three new options appear here: `"Device" -> "TBB"` (CPU rather than GPU),
`"AD" -> True` (emit automatic-differentiation instantiations too), and empty `"Coordinates"`
(no external-momentum grid — the potential is a function of the field, not of $p$).

The resulting `Vsig` kernel is the smallest real kernel in the whole tree — 518 bytes, `nenv = 0`,
one trace function that returns the literal `1` — and its body is the Wetterich equation in a line:

```cpp
const auto _den1 = powr<-1>(_interp1 + d1VL + _interp3*powr<2>(l1) + 2.*d2VL*rhoL);
return 0.5 * fma(_den1, tr0(fenv)*_interp2, fma(-1., _den1*tr0(fenv)*_interp1*etaSigL, 0.));
```

That is $\tfrac12\int (\dot R - \eta R)/(Z p^2 + m^2 + R)$. Worth looking at once: it is the whole
formalism, emitted.

## Possibilities for extensions

1. **Derive `spAngles3`.** From the `sp3Frame` table in [step-07](step-07.md), compute
   $\hat l\cdot\hat q_2$ and confirm the expression above. Then derive `spAngles4` from `sp4Frame`
   and check it against the production file.

2. **Break the frame/integrator pairing.** Change `"Integrator"` to `Integrator_p2_1ang` while
   keeping `frame2`. Find out how far the mistake gets before anything complains.

3. **Watch a restriction pay.** Generate `ZAqbq1` with and without `dropAqbq47` and compare diagram
   counts, generation time and emitted size. Then reproduce the 133/439/1349 numbers by setting
   `NT_FQ_TENSORS`.

4. **Cut in the wrong place.** Write `dropAqbq47` (named heads) but apply it *before*
   `dressingRules`. The rule never fires; the kernel is the unrestricted one. Confirm that the only
   symptom is a size difference.

5. **Compare the two symmetry spellings.** Run `ZA3` with `FMakeSymmetryList[{A[i1],A[i2],A[i3]}]`
   and with `FSymmetry[Symmetric, {i1,i2}]` and compare diagram counts *and values*. They must agree
   numerically; if they do not, the full group was not a symmetry.

6. **Reconstruct the potential.** Verify `V == Vadv + Vsig` numerically from the three emitted
   kernels. It should hold identically, since the supertrace slices are disjoint.

Next, [step-15](step-15.md): what `MakeNTKernelDiFfRG` actually wrote to disk, and the contract
between that `.wl` parameter list and the C++ that calls it.
