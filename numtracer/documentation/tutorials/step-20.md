# step-20: Making a hard flow generate

*Builds on: [step-09](step-09.md), [step-15](step-15.md) · Built on by: [step-22](step-22.md) · Tags: `performance`, `build` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* Phase A (parallel Dirac contraction, RAM ∝ workers) versus phase B (per-net fold *and* the CSE/Horner lowering) — and why throttling phase B globally is a mistake.
* `SetNumTracerThreads[nA, nB]`, `NT_GEN_MAXW` / `NT_GEN_MAXW_B`, and a RAM budget.
* Chunking (`ntChunkDef`) so a single emitted function never blows the compiler; `NT_GEN_MAIN_OPT` and why `-O2` on the main TU can be 100x worse than `-O1`.
* Sub-term dedup across diagrams, and `"Offline" -> True` + `numtracer_add_numtrace` so generation happens in the consumer's build.
* The full `NT_GEN_*` knob catalogue, and which ones are measured losses.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `numtracer/PERFORMANCE.md`
* the root `NUMTRACER_*.md` design briefs
