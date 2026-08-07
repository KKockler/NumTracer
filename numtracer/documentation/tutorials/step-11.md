# step-11: Projectors and tensor bases

*Builds on: [step-08](step-08.md) · Built on by: [step-12](step-12.md), [step-14](step-14.md) · Tags: `tensorbases`, `projector` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* `TBGetProjector` vs `TBGetVertex` — a projector always projects *vertices*; never reach for `TBGetInnerProduct`.
* Index-order and sign conventions, and the momentum-routing change that once reversed the quark loop.
* Asserting a projector instead of trusting it: the `must be 1 … else Abort[]` self-check pattern from `QCD_Nf2/with_mesons`.
* Orthonormality as a physics gate: `<P_i, T_j> = delta_ij` over a whole basis, needing no tracer at all.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `numtracer/tests/gen/gen_ortho_numeric.wls` and `tests/compare_ortho_num.cpp`
* `QCD_Nf2/with_mesons/QCD.wl`, the `ProjectorhPi` block
