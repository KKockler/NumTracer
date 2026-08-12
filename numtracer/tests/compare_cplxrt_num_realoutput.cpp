// The "RealOutput" half of the finite-density gate (see compare_cplxrt_num.cpp).
//
// gen_cplxrt_numeric.wls emits the SAME traced kernel twice: once normally (three bodies) and once
// with "RealOutput" -> True, which drops the untouched complex body a real-valued consumer could
// never instantiate. This TU compiles the RealOutput kernel with the verdict forced to 0 — the case
// that option exists to define. With no complex body to fall through to, verdict 0 selects the RePart
// body behind a #warning, i.e. the kernel deliberately returns Re[...] of a genuinely complex flow.
//
// That is a TRUNCATION of the flow equation, not an identity, and the point of grading it here is
// that the truncation must still be the exact real part: whatever the physics of dropping Im, the
// number returned has to equal real(complex body) pointwise. If it does not, the RealOutput branch is
// not selecting the body it claims to.
//
// The expected #warning is suppressed for this TU only (see CMakeLists.txt): it is the whole point of
// the file, and a warning fired on purpose should not look like a defect in the build log.
#include "numtrace_verdict.hh"
#undef NT_CPLXRTRO_NUM_VERDICT
#define NT_CPLXRTRO_NUM_VERDICT 0

#include "Cplxrtro_num_kernel.hh"

double cplxrt_realoutput(double l0, double l1, double cos1, double p0, double p, double muq,
                         double Ep)
{
  return numtracer_kernels::Cplxrtro_num_kernel::kernel(l0, l1, cos1, p0, p, muq, Ep);
}
