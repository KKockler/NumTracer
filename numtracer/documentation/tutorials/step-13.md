# step-13: Counterterms and the constant slot

*Builds on: [step-12](step-12.md) · Built on by: [step-14](step-14.md) · Tags: `frg`, `diffrg` ·
**Tier C** (the full fRG toolchain)*

```{admonition} What this step needs
:class: important
The material is read from `QCD_Nf2/baseline/QCD.wl` and `QCD_Nf2/no_mesons/model.hh`. Nothing here
is a new script — this step is about a mechanism you have already seen emitted
([step-09](step-09.md)) and about the C++ that consumes it.
```

## Introduction

[step-09](step-09.md) showed that `"Constant" -> expr` emits a second entry point:

```cpp
static inline auto constant(const double& p, const auto& myZ) { return myZ(p); }
```

and that it leaves `kernel()` untouched. That was the mechanism. This step is *why anyone would
want it* — and the answer turns out to be load-bearing: without it, the quark sector of a QCD flow
does not converge at all.

### Anomalous dimensions and the pinning condition

A propagator dressing $Z_q(p)$ is defined only up to an overall constant: rescaling the quark field
rescales $Z_q$ everywhere. Numerically you must therefore **fix a normalisation**, and the usual
choice is to pin the dressing at one momentum,

$$
Z_q(p_{\min}) \;=\; 1 \quad\text{for all RG time } t .
$$

Holding that fixed as $k$ runs means continuously rescaling the field, and the rate of rescaling is
the **anomalous dimension** $\eta_q$. It enters the flow as an extra term:

$$
\partial_t Z_q(p) \;=\; \underbrace{[\text{loops}]}_{\text{the integral}} \;+\; \eta_q\,Z_q(p).
$$

That second term is the counterterm. Two things about it:

* **It is loop-independent.** $\eta_q Z_q(p)$ does not depend on the loop momentum, so integrating
  it over the loop would be wasted work — worse, it would be multiplied by the measure.
* **It is what makes the fixed point solvable.** $\eta_q$ is not known in advance; it is
  *determined* by demanding $\partial_t Z_q(p_{\min}) = 0$. That is one scalar equation in one
  unknown, solved by a Newton iteration at every step — and the counterterm is the term that
  supplies the diagonal.

Hence a separate entry point: something the integrator adds *after* quadrature.

## The Mathematica side

From `QCD_Nf2/baseline/QCD.wl`, the `Zq` flow:

```mathematica
tr  = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRG["1-Loop"]["ExternalIndices"]],
            (fRG["1-Loop"]["Expression"] /. diagRules)];
ntk = NumTrace[FromFunKit[PropParam[dressingRules[tr]]], "Frame" -> frame1, "Args" -> args1];
MakeNTKernelDiFfRG[ntk, etaQ*zq[p],          (* <-- positional constExpr, 2nd argument *)
    "Name" -> "Zq", "Namespace" -> "zq_qcd",
    "Integrator" -> "Integrator_p2_1ang",
    "Parameters" -> kernelParameterList,
    "IntegrationVariables" -> {"l1", "cos1"}];
```

`etaQ*zq[p]` sits in the **positional second argument**. It is a *plain expression*, not an
`NTKernel` — it never goes near the tracer. `MakeNTKernelDiFfRG` accepts it either positionally as
above or as `"Constant" -> expr`.

The production comment states the contract:

```{admonition} From QCD_Nf2/baseline/QCD.wl
:class: note
> The positional `constExpr` `etaQ*zq[p]` is the RG counterterm that pins `zq(p_min) = 1`; it
> becomes `constant()` in the emitted kernel and supplies the **+1 diagonal the Newton loop needs**.
> `ZAqbq1` carries it too.
```

Note that `etaQ` is a *scalar kernel parameter*, declared in `$ntScalars` alongside `k`:

```mathematica
$ntScalars = {
   <|"Name" -> "k",    "Type" -> "double"|>,
   <|"Name" -> "etaQ", "Type" -> "double"|>};
```

It is a number the C++ passes in and updates between iterations — not a dressing, not an
interpolator.

### The four-fermi generalisation

A counterterm's coefficient counts the field rescalings the object carries. A four-quark coupling
has *two* quark bilinears, so $\bar\lambda = \lambda/Z_q^2$ and the counterterm doubles:

```mathematica
MakeNTKernelDiFfRG[ntk, 2*etaQ*Symbol["lambda4L" <> ToString[n]][p],
    "Name" -> "lambda4L" <> ToString[n], …];
```

And a hadronised Yukawa carries one quark bilinear plus one meson leg, so it gets one of each —
plus the hadronisation term itself:

```mathematica
MakeNTKernelDiFfRG[ntk, (etaPiL/2 + etaQ)*hPhiL[S0, S1, SPhi] + hadronisationA6, …]
```

The `constExpr` slot takes an arbitrary expression, not just a single product.

## The C++ side

The emitted kernel gains its second entry point:

```cpp
static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const double& etaQ,
                                            const SplineInterpolator1D<…>& zq, …)
{
  const auto _interp1 = zq(p);
  return _interp1 * etaQ;
}
```

`constant()` takes the coordinate arguments and the parameters, but **not** the integration
variables — there is nothing to integrate. DiFfRG's integrator calls it once per external momentum
and flat-adds the result to the quadrature.

### The Newton loop that consumes it

From `QCD_Nf2/no_mesons/model.hh` — the self-consistent solve for the anomalous dimensions, run at
every RG step *before* the actual residuals are computed:

```cpp
bool eta_converged = false;
int n_iter = 0;
while (!eta_converged) {
  for (uint i = 0; i < p_grid_size; ++i) { old_dtZA[i] = dtZA[i]; old_dtZc[i] = dtZc[i]; }
  flow_equations.ZA.map(&residual[idxv("ZA")], coordinates1D, arguments);
  flow_equations.Zc.map(&residual[idxv("Zc")], coordinates1D, arguments);
  dtZA.update(&residual[idxv("ZA")]);
  dtZc.update(&residual[idxv("Zc")]);

  double dist = 0.;
  for (uint i = 0; i < p_grid_size; ++i) {
    dist = std::max(dist, std::abs(dtZA[i] - old_dtZA[i]) / std::abs(dtZA[i]));
    dist = std::max(dist, std::abs(dtZc[i] - old_dtZc[i]) / std::abs(dtZc[i]));
  }

  // quark: pin zq(p_min) = 1 via etaQ (Newton on the undamped residual, under-relaxed step)
  flow_equations.Zq.map(zqF.data(), coordinates1D, arguments);
  const double resQ = zqF[0] / variables[idxv("zq")];
  etaQ -= prm.eta_relax * resQ;
  dist = std::max(dist, std::abs(resQ) / (std::abs(etaQ) + 1e-30));

  n_iter++;
  if (dist < prm.eta_tol || n_iter >= prm.eta_iter_max) eta_converged = true;
}

// quark dressing: flows with the +etaQ*zq counterterm baked into its kernel
flow_equations.Zq.map(&residual[idxv("zq")], coordinates1D, arguments);
```

Read `resQ` carefully. `zqF[0]` is $\partial_t Z_q$ at the **lowest grid momentum** — including the
counterterm, because it is baked into the kernel. Dividing by $Z_q$ there and demanding it vanish
*is* the pinning condition. The update `etaQ -= eta_relax * resQ` is a Newton step, under-relaxed
because the map is stiff.

**Why the counterterm makes this converge.** Without it, $\partial_t Z_q(p_{\min})$ would be the
loop integral alone, with no explicit dependence on $\eta_q$ — the Newton iteration would have no
derivative to work with. With it, $\partial_t Z_q = [\text{loops}] + \eta_q Z_q$, so
$\partial(\text{resQ})/\partial\eta_q \approx 1$: the "+1 diagonal" the comment refers to. The
iteration is then a contraction near the solution.

Note also the loop is a *coupled* solve — $\eta_A$, $\eta_c$ and $\eta_q$ are found together,
because each flow's regulator insertion depends on the others' $\partial_t Z$.

### When it goes wrong

This is why [step-12](step-12.md)'s normalisation rule matters so much:

```{admonition} From QCD_Nf2/baseline/QCD.wl
:class: warning
> an extra `1/sp[p1, p1]` would over-normalise the flow by $1/p^2$. At `p_grid_min = 1e-3` that is a
> **1e6 enhancement** at the pinning point, which turns the `etaQ` Newton iteration in `model.hh`
> into an **expanding map** (`etaQ -> 1e250 -> NaN` on the first timestep).
```

The counterterm is $+\eta_q Z_q$, i.e. a derivative of 1. If the loop part is over-normalised by
$10^6$ *at the pinning momentum*, that 1 is swamped, the Newton step overshoots by six orders of
magnitude, and the iteration diverges immediately. A projection prefactor and a Newton iteration
that look completely unrelated are in fact the same equation.

## Possibilities for extensions

1. **Find the diagonal.** Emit the `Zq` kernel with and without the `constExpr`, and evaluate
   `constant()` and `kernel()` separately at $p_{\min}$. Confirm numerically that
   $\partial(\text{resQ})/\partial\eta_q \approx 1$ with the counterterm and $\approx 0$ without.

2. **Break the relaxation.** Set `eta_relax` to 1.0 (an undamped Newton step) and watch where it
   stops converging. Then reproduce the $10^6$ disaster by adding a spurious `1/sp[p1,p1]` to the
   `Zq` projection, and observe that no amount of relaxation saves it.

3. **Count the rescalings.** For each flow in `with_mesons` that carries a `constExpr`, check that
   the coefficient of `etaQ` equals the number of quark bilinears in the object. `lambda4L*` has 2;
   `ZAqbq1` has 1; `hPhiL` has 1 plus a meson factor. A flow with the wrong count has a wrong
   anomalous dimension and will drift.

4. **Look at a `constant()` that is zero.** Most flows do not carry one, and
   `MakeNTKernelDiFfRG` emits `return 0.;`. Confirm that DiFfRG still calls it, and think about why
   emitting it unconditionally is better than making the interface optional.

5. **Trace the parameter.** Follow `etaQ` from `$ntScalars` in the `.wl`, through
   `kernelParameterList`, into the emitted signature, into `device::tie(k, etaQ, …)` in `model.hh`.
   Then reorder `$ntScalars` and see how far the mistake travels before anything complains —
   [step-15](step-15.md).

Next, [step-14](step-14.md): three- and four-point vertices, where the frame grows angles and one
derivative feeds several projectors.
