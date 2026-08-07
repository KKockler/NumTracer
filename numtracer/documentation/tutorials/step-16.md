# step-16: Finite temperature

*Builds on: [step-04](step-04.md), [step-07](step-07.md) · Built on by: — · Tags: `finite-T`, `projector` · **Tier A + B***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* The heat bath breaks O(4) to spatial O(3): component 0 becomes physically distinct.
* `ntElectricProj` / `ntMagneticProj` — the thermal split of $P^T$, with $\mathrm{tr}P^E = 1$ and $\mathrm{tr}P^M = 2$.
* `ntSPS` (spatial scalar product) and `ntVec[q, 0]` (an *integer* index = the Matsubara component).
* `propFrameFT`, and why an independent temporal loop component takes the general frame path.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `numtracer/tests/test_finiteT_quark.cpp`
* `numtracer/tests/gen/gen_ftproj_numeric.wls`
