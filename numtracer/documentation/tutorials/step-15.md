# step-15: Wiring a flow into DiFfRG

*Builds on: [step-14](step-14.md) · Built on by: [step-20](step-20.md), [step-22](step-22.md) · Tags: `diffrg`, `build` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* `MakeNTKernelDiFfRG` and `UpdateNTFlows`, and why the two are bundled.
* The anatomy of a generated `flows/<Name>/`: `kernel.hh`, `kernels.hh`, `<Name>.hh`, `src/`, `numtrace.json`, `sources.m`.
* `numtracer_add_numtrace` and the offline `numtrace` build target.
* **The parameter-order contract**: `kernelParameterList` in the `.wl` *is* the C++ call signature, and nothing checks it.
* How `model.hh` consumes a kernel: `device::tie(...)`, `.map(...)`, and the interpolator update dance.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `QCD_Nf2/baseline/` in full — 588-line `QCD.wl`, 245-line `model.hh`, `flows/CMakeLists.txt`
