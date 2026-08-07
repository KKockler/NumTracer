# step-22: Validating and debugging a kernel

*Builds on: [step-15](step-15.md), [step-20](step-20.md) · Built on by: — · Tags: `testing`, `debugging` · **Tier C***

```{admonition} This page is not written yet
:class: warning
The tutorial series is being written in order. This step is planned and its scope is
fixed, but the narrative below is an outline rather than a finished page.
```

## What this step will cover

* **Plain `ctest` is blind to codegen.** It compiles pre-generated headers and never invokes Wolfram, so you can break the generator completely and still see all tests pass. `-DNUMTRACER_TEST_CODEGEN=ON` and `ctest -L codegen` is the real gate.
* Oracles: what `tests/gen/` freezes, and the rule that both sides must be frozen from the same run.
* Identity checks that need no oracle: projector orthonormality, $P^T + P^L = \delta$, a general frame reducing to a special one.
* Reading a failure — the signatures worth recognising: a ratio of exactly 0.25 (a collapsed Dirac loop), a leaked `dressing(...)` call, a `$`-suffixed symbol in generated code, a stale archive giving an ODR mismatch.

## Where the material lives in the meantime

Until this page is written, read these directly:

* `numtracer/tests/gen/README.md`
* `STALE_TEST_ORACLES.md`
