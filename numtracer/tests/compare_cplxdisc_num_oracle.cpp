// The ORACLE half of the complex x multi-trace gate (see compare_cplxdisc_num.cpp).
//
// MakeNTKernel emits all three bodies into one header behind `#if NT_CPLXDISC_NUM_VERDICT`, so only
// one of them is compiled per translation unit. This TU compiles the UNTOUCHED COMPLEX body — the one
// that performs every complex multiplication in C++ and is therefore correct by construction, at any
// token degree — and hands its real part to the main TU, which compiles the RePart body the committed
// verdict selects. Comparing the two is the whole test.
//
// The macro must be forced BEFORE the kernel header pulls in the generated numtrace_verdict.hh.
// Including that header first is what makes this legal rather than a redefinition: it is `#pragma
// once`, so the kernel header's own `#include "numtrace_verdict.hh"` is then a no-op and the #undef
// below is not overwritten.
#include "numtrace_verdict.hh"
#undef NT_CPLXDISC_NUM_VERDICT
#define NT_CPLXDISC_NUM_VERDICT 0

#include "Cplxdisc_num_kernel.hh"

#include <complex>

// A distinct symbol, not the kernel class itself: the main TU compiles a DIFFERENT body under the
// same class name, so exposing the class from here would be an ODR violation. Only this plain
// function crosses the TU boundary.
double cplxdisc_oracle(double l1, double cos1, double cos2, double p)
{
  return std::real(numtracer_kernels::Cplxdisc_num_kernel::kernel(l1, cos1, cos2, p));
}
