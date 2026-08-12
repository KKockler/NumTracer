// The ORACLE half of the finite-density gate (see compare_cplxrt_num.cpp).
//
// MakeNTKernel emits all three bodies into one header behind `#if NT_CPLXRT_NUM_VERDICT`, so only one
// of them is compiled per translation unit. This TU compiles the UNTOUCHED COMPLEX body — the one that
// performs every complex multiplication in C++, including the complex denominator, and is therefore
// correct by construction. Its real part is the oracle for the two real bodies compiled elsewhere.
//
// The macro must be forced BEFORE the kernel header pulls in the generated numtrace_verdict.hh.
// Including that header first is what makes this legal rather than a redefinition: it is `#pragma
// once`, so the kernel header's own `#include "numtrace_verdict.hh"` is then a no-op and the #undef
// below is not overwritten.
#include "numtrace_verdict.hh"
#undef NT_CPLXRT_NUM_VERDICT
#define NT_CPLXRT_NUM_VERDICT 0

#include "Cplxrt_num_kernel.hh"

#include <complex>

// A distinct symbol, not the kernel class itself: the other TUs compile DIFFERENT bodies under the
// same class name, so exposing the class from here would be an ODR violation. Only this plain
// function crosses the TU boundary.
double cplxrt_oracle(double l0, double l1, double cos1, double p0, double p, double muq, double Ep)
{
  return std::real(
      numtracer_kernels::Cplxrt_num_kernel::kernel(l0, l1, cos1, p0, p, muq, Ep));
}
