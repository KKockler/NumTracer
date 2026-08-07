# step-12: A propagator flow, end to end

*Builds on: [step-10](step-10.md), [step-11](step-11.md) · Built on by: [step-13](step-13.md) · Tags: `frg`, `flow` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* The ghost and gluon two-point flows, `Zc` and `ZA`, from derivative to kernel.
* Leg symmetries: `FMakeSymmetryList[FSymmetry[Symmetric, {i1,i2}]]` versus the field-list form, and *when the reduction inside `FTruncate` is legitimate*.
* Reading the four terms of the emitted `ZA` kernel as four physical diagrams.
* The compensating `-1/sp[p1,p1]` on the ghost projection, and what dropping it does.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `QCD_Nf2/baseline/QCD.wl` and `baseline/flows/{Zc,ZA}/`
