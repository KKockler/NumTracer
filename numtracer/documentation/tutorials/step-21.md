# step-21: GPU kernels and complex traces

*Builds on: [step-09](step-09.md) · Built on by: — · Tags: `gpu`, `cuda`, `complex` · **Tier B***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* `"Decorator" -> "static __host__ __device__ inline"` — the same straight-line arithmetic, now a device function.
* The two-phase quadrature integrator and its Kokkos twin.
* Why a trace is complex at all, and the three emitted kernel bodies (untouched / pure-imaginary-dropped / real-part).
* The generated probe that picks one and writes `numtrace_verdict.hh` — and why `PruneRealTraces` is hazardous in combination with it.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `numtracer/tests/gpu/README.md` and `tests/gpu/integrator.cuh`
