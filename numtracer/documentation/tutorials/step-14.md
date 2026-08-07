# step-14: Three- and four-point vertices

*Builds on: [step-11](step-11.md), [step-13](step-13.md) · Built on by: [step-15](step-15.md), [step-19](step-19.md) · Tags: `vertex`, `frame` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* `ZAcbc` -> `ZA3` -> `ZAqbq1/4/7`: three distinct fields, then full Bose symmetry, then one derivative with several projectors.
* `SPParam` and `"AngleDefs" -> spAngles3` — keeping the scalar and tensor halves on the same kinematics ([step-07](step-07.md)).
* `TBRestrictBasis`, and memoising a shared derivative (`zaqbqFRG := zaqbqFRG = …`).
* Restricting a trace with `dropAqbq47` / `dropFourQuark`, and what it costs.
* As an aside: slicing the supertrace directly (`(WetterichEquation//FTruncate)[[{3,4}]]`) for an effective-potential flow.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `QCD_Nf2/baseline/QCD.wl`, `no_mesons/QCD.wl`
