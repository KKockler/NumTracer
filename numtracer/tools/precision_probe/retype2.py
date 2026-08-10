#!/usr/bin/env python3
"""G1 retype, v2 — supersedes retype.py (which predates the `nt_complex_t`/`NT_TRACE_COMPLEX`
element type introduced in 954de3e and no longer matches the emitted kernels).

Two scopes:

  traces : the emitted SSA layer only (kernels.hh powr/fill/trN + the fenv array). Splines,
           `_interp`/`_den` defs and the assembly stay fp64. This is what retype.py did.
  full   : additionally retypes the wrapper's arithmetic — cosl1p*, dr_*, and the _interp/_den
           definitions — so the whole per-node computation runs in the reduced type. Spline
           STORAGE stays fp64 (it is DiFfRG-side data, and its return value converts on
           assignment); the kernel signature also stays double, since the integrator calls it
           with doubles and nt_real converts implicitly.

usage: retype2.py <src_dir(tests/gen)> <dst_dir> [traces|full]
"""
import re
import sys
from pathlib import Path

src_dir, dst_dir = Path(sys.argv[1]), Path(sys.argv[2])
scope = sys.argv[3] if len(sys.argv) > 3 else "traces"
assert scope in ("traces", "full"), scope
dst_dir.mkdir(parents=True, exist_ok=True)


def _d(expr):
    """Wrap a returned SSA operand in the explicit nt_real -> double conversion, unless it is
    already a plain literal (0.0 / -0.0), which converts on its own."""
    e = expr.strip()
    return e if re.fullmatch(r"-?\d+(\.\d*)?", e) else "nt_to_double(%s)" % e

# ---- traces header -------------------------------------------------------------------------
ks = (src_dir / "ZA3_gpu_kernels.hh").read_text()
ks = ks.replace("#include <cmath>", '#include <cmath>\n#include "nt_real.cuh"', 1)
# Drop the emitted `#ifndef NT_TRACE_COMPLEX / #define ... std::complex<double>` default: the GPU
# shim already defines NT_TRACE_COMPLEX as cuda::std::complex<double>, and leaving the default here
# would expose a `double` to the rewrite below.
ks = re.sub(r"#ifndef NT_TRACE_COMPLEX\n#define NT_TRACE_COMPLEX [^\n]*\n#endif\n", "", ks)
# Retype the SSA INTERNALS only. The trace INTERFACE (return types, and nt_complex_t itself) stays
# as emitted, because tests/gpu/shim.hpp deliberately unifies the traces' complex type and the COEN
# assembly's complex coefficients on cuda::std::complex<double> — retyping the return type splits
# them again and there is no operator+ across the two (shim.hpp:52-64).
ks = re.sub(r"\bdouble\b", "nt_real", ks)
# ...so restore the emitted signatures, and convert at the return boundary (nt_real's operator
# double is explicit by design, so the conversion has to be written out).
ks = re.sub(r"(KOKKOS_[A-Z_]*FUNCTION )nt_real (tr\d+|powr)", r"\1double \2", ks)
ks = re.sub(r"\breturn ([A-Za-z_]\w*);", r"return nt_to_double(\1);", ks)
ks = re.sub(r"return nt_complex_t\{([^,}]+), ([^}]+)\};",
            lambda m: "return nt_complex_t{%s, %s};" % (_d(m.group(1)), _d(m.group(2))), ks)
(dst_dir / "ZA3_gpu_kernels.hh").write_text(ks)

# ---- kernel wrapper ------------------------------------------------------------------------
k = (src_dir / "ZA3_gpu_kernel.hh").read_text()
k = k.replace("double fenv[", "nt_real fenv[")
if scope == "full":
    # body scalars only — never the signature (spline params must stay SplineInterpolator1D<double>)
    k = re.sub(r"\bconst double (cosl1p\d+|dr_\d+)\b", r"const nt_real \1", k)
    # force the hoisted interpolator/denominator temporaries into the reduced type; the spline
    # call still returns double and converts on assignment
    k = re.sub(r"\bconst auto (_interp\d+|_den\d+|_cse\d+)\b", r"const nt_real \1", k)
(dst_dir / "ZA3_gpu_kernel.hh").write_text(k)

cu = (src_dir.parent / "gpu" / "integrate_za3_gpu.cu").read_text()
(dst_dir / "integrate_za3_gpu_variant.cu").write_text(cu)
print(f"retyped ({scope}) into {dst_dir}")
