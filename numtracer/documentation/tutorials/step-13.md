# step-13: Counterterms and the constant slot

*Builds on: [step-12](step-12.md) · Built on by: [step-14](step-14.md) · Tags: `frg`, `diffrg` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* The RG counterterm that pins a dressing at the lowest grid point, and why `Zq` carries *no* `1/sp` prefactor.
* `MakeNTKernelDiFfRG[ntk, etaQ*zq[p], …]` — the positional `constExpr` slot, and the `constant()` it becomes ([step-09](step-09.md) showed the mechanism).
* The Newton iteration in `model.hh` that consumes it, and the 1e6-enhancement-to-NaN failure mode when the normalisation is wrong.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `QCD_Nf2/baseline/QCD.wl` (Zq) and `baseline/model.hh`
