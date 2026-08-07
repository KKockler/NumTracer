# step-18: Dressed numerators and vertex collection

*Builds on: [step-03](step-03.md), [step-05](step-05.md) · Built on by: — · Tags: `dressing`, `dpoly` · **Tier A + B***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* $S(p) = M_q\mathbb{1} + Z(p)\slashed p$ kept eager: one collected trace instead of $2^D$ diagrams.
* `DSlot` / `DSlotOpt` / `DChainTok`, and the `DPoly` the collection folds into.
* The general `toks` / `netFacs` form, which also covers open-leg *vertex* structures.
* Vertex collection (`ntDiracSlot`, `NT_VERTEX_COLLECT=1`) and when it wins.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `Tutorials/step-18-dressed-numerators/dressed_numerators.cpp` (runnable now: `ctest -R dressed_numerators`)
