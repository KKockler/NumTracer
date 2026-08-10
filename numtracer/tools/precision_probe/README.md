# G1 precision-ladder probe (2026-08-08)

Measures GPU runtime/accuracy of the emitted trace layer in reduced precision, exploiting
sm_89's 1:64 fp64:fp32 throughput ratio. `nt_real.cuh` provides the wrapped scalar
(NT_REAL_MODE=1 fp32, =2 double-single) — a wrapper because the kernels are full of double
literals that would silently re-promote bare float arithmetic. `retype.py` rewrites copies of
ZA3_gpu_kernel(s).hh onto it (traces+env only; splines/assembly stay fp64).

## 2026-08-08 update — gate is GREEN, re-measured on the regenerated (954de3e) kernels

`retype.py` is SUPERSEDED by `retype2.py`: the emitted traces header now types its complex returns
via `nt_complex_t = NT_TRACE_COMPLEX`, so the old `std::complex<double>` substitution no longer
matches anything. `retype2.py` also keeps the trace **interface** (return types, `nt_complex_t`) as
emitted and retypes only the SSA internals — `tests/gpu/shim.hpp:52-64` deliberately unifies the
traces' complex type with the COEN assembly's complex coefficients on `cuda::std::complex<double>`,
and splitting them leaves no `operator+` between the two.

The blocker is gone: `ref_za3` now reads **1.003e-15 PASS** (both sides regenerated together).

| variant | phase1 | speedup | max rel err |
|---|---|---|---|
| fp64 baseline | 3.956 ms | 1.00x | 1.003e-15 |
| double-single traces (`NT_REAL_MODE=2`) | 3.309 ms | **1.20x** | 5.732e-08 |
| fp32 traces (`NT_REAL_MODE=1`) | 2.738 ms | **1.44x** | 8.904e-06 |

Both miss the harness's 1e-8 gate — double-single only just (5.7e-08). That gate was chosen to
separate a mechanics BUG (<1e-8) from legitimate FMA/summation reordering (~1e-13), so 5.7e-08 is a
real precision loss, not a defect; whether it is acceptable is a physics call about the flow's
error budget, not a code question.

Note the ratios IMPROVED vs the pre-regeneration numbers below (1.28x/1.07x): the regenerated
kernels put more work in the trace layer. The "full retype approaches 2x" hypothesis now looks
unreachable — traces-only already captures 1.44x, and the remainder is spline lookup and the
fp64 assembly. `retype2.py --scope full` exists but is UNTESTED.

### Original (pre-954de3e) numbers, kept for the record

Build (from numtracer/):
  python3 tools/precision_probe/retype.py tests/gen /tmp/v32 && cp tools/precision_probe/nt_real.cuh /tmp/v32/
  nvcc -arch=sm_89 -O3 --expt-relaxed-constexpr -DNT_REAL_MODE=1 -I /tmp/v32 -I tests/gpu \
       -I tests/gen -I include /tmp/v32/integrate_za3_gpu_variant.cu build/libnt_gpu_ref.a -lgsl -lgslcblas -o /tmp/za3_fp32
