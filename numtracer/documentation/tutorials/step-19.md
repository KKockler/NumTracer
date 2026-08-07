# step-19: Four-fermi, Fierz bases, fixed-component gammas

*Builds on: [step-08](step-08.md), [step-14](step-14.md) · Built on by: — · Tags: `four-fermi`, `fierz` · **Tier B***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* Two *independent* spinor lines in one diagram (`dloopsep`), and momentum-independent contact structures.
* Fixed-component gammas: $\gamma^0$ is a slash against the unit basis vector $e_0$ (`ntUnitVec`), needing no new engine machinery.
* The Fierz-complete four-quark basis and its orthonormality gate.
* The symmetry subtlety: `FSymmetry[Symmetric,{i1,i3},{i2,i4}]` is respected by 78/78 diagrams where the auto field-list form is respected by only 14/78.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `numtracer/tests/gen/gen_fierz_ortho_numeric.wls`, `gen_fixedcomp_numeric.wls`
