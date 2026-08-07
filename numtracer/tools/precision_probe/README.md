# G1 precision-ladder probe (2026-08-08)

Measures GPU runtime/accuracy of the emitted trace layer in reduced precision, exploiting
sm_89's 1:64 fp64:fp32 throughput ratio. `nt_real.cuh` provides the wrapped scalar
(NT_REAL_MODE=1 fp32, =2 double-single) — a wrapper because the kernels are full of double
literals that would silently re-promote bare float arithmetic. `retype.py` rewrites copies of
ZA3_gpu_kernel(s).hh onto it (traces+env only; splines/assembly stay fp64).

Measured on ZA3_gpu (best of 3, phase1): fp64 1.59 ms, fp32 traces 1.24 ms (1.28x),
double-single 1.48 ms (1.07x); integral samples agree to ~7 digits (fp32 eps level).
A full-kernel retype could approach 2x, but the accuracy gate needs a REFRESHED long-double
reference first — ref_za3 is stale at 3.6e-02 even against the committed fp64 kernels.

Build (from numtracer/):
  python3 tools/precision_probe/retype.py tests/gen /tmp/v32 && cp tools/precision_probe/nt_real.cuh /tmp/v32/
  nvcc -arch=sm_89 -O3 --expt-relaxed-constexpr -DNT_REAL_MODE=1 -I /tmp/v32 -I tests/gpu \
       -I tests/gen -I include /tmp/v32/integrate_za3_gpu_variant.cu build/libnt_gpu_ref.a -lgsl -lgslcblas -o /tmp/za3_fp32
