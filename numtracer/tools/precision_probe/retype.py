#!/usr/bin/env python3
"""G1: retype the ZA3_gpu trace layer (kernels.hh: powr/fill/trN SSA) to nt_real, and the
fenv array in kernel.hh to match. Spline calls, _interp defs and the assembly stay double —
the probe measures the trace arithmetic, which is the emitted-SSA bulk."""
import re
import sys
from pathlib import Path

src_dir, dst_dir = Path(sys.argv[1]), Path(sys.argv[2])
dst_dir.mkdir(parents=True, exist_ok=True)

ks = (src_dir / "ZA3_gpu_kernels.hh").read_text()
ks = ks.replace("#include <complex>", "#include <complex>\n#include \"nt_real.cuh\"")
ks = ks.replace("std::complex<double>", "nt_complex")
ks = re.sub(r"\bdouble\b", "nt_real", ks)
# the env-size constant must stay integral
ks = ks.replace("static inline constexpr int nenv", "static inline constexpr int nenv")
(dst_dir / "ZA3_gpu_kernels.hh").write_text(ks)

k = (src_dir / "ZA3_gpu_kernel.hh").read_text()
k = k.replace("double fenv[", "nt_real fenv[")
(dst_dir / "ZA3_gpu_kernel.hh").write_text(k)

# harness: same TU, just picks up the variant headers via -I order
cu = (src_dir.parent / "gpu" / "integrate_za3_gpu.cu").read_text()
(dst_dir / "integrate_za3_gpu_variant.cu").write_text(cu)
print("retyped into", dst_dir)
