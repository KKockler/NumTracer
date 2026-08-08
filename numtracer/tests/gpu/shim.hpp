// NumTracer — host/device shim for the GPU integration tests. The SAME header serves
//   * the .cu TU evaluating the __host__ __device__ generated kernels on the GPU, and
//   * the plain C++ TU computing the CPU same-quadrature reference with the host kernels
// (the GPU test targets put tests/gpu BEFORE tests/refshim on the include path, so the
// generated kernels' `#include "shim.hpp"` resolves here for both TUs).
//
// Every regulator/dressing FORMULA is copied verbatim from tests/refshim/shim.hpp — the
// GPU validation compares against the CPU path, so both sides must be the same algebraic
// functions. Structural differences only:
//   * NT_DH expands to __host__ __device__ under a CUDA compiler, else to nothing;
//   * Fn dispatches on a runtime id instead of wrapping a host function pointer (a host
//     fn-pointer cannot be called on the device, and all 11 dressings share the single
//     SplineInterpolator1D alias type, so tag-type dispatch is impossible).
#pragma once

#include <cmath>
#include <complex>

#if defined(__CUDACC__) || defined(__CUDA__)
#define NT_DH __host__ __device__
#include <cuda/std/complex>
#else
#define NT_DH
#endif

#define KOKKOS_FORCEINLINE_FUNCTION inline
// The generated GPU kernels now decorate with KOKKOS_INLINE_FUNCTION (the gpu .wls emit for the
// Kokkos consumer); this TU is not a Kokkos build, so supply the equivalent expansion. Guarded:
// in a real Kokkos TU the macro is already defined and must win.
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION NT_DH inline
#endif

// ---- the regulators the generated kernels call UNQUALIFIED -------------------
// NumTracer emits a plain kernel class; the consumer supplies RB/RF/... . Global scope,
// so unqualified lookup from a kernel member (class -> namespace -> global) finds them.
// Same formulas as tests/refshim/nt_regulators.hpp, carrying NT_DH because the GPU
// kernels are __host__ __device__ and a plain host function cannot be called there.
NT_DH inline double RB(double k2, double p2) { return k2 * std::exp(-p2 / k2); }
NT_DH inline double RBdot(double k2, double p2) { return 2.0 * (k2 + p2) * std::exp(-p2 / k2); }
NT_DH inline double RF(double k2, double p2) { return std::sqrt(k2) * std::exp(-p2 / k2); }
NT_DH inline double RFdot(double k2, double p2) { return (k2 + 2.0 * p2) / std::sqrt(k2) * std::exp(-p2 / k2); }
NT_DH inline double dq2RB(double, double) { return 0.0; }
NT_DH inline double dq2RF(double, double) { return 0.0; }

namespace DiFfRG {

// The complex type the generated kernels resolve UNQUALIFIED (they are emitted without
// `using std::complex;`): std::complex's arithmetic lowers to gcc `_Complex` builtins that
// nvcc silently miscompiles to 0 in device code, so CUDA TUs must use libcu++'s
// cuda::std::complex instead. Host TUs keep std::complex.
// RE-CONFIRMED 2026-08-08 on CUDA 12.9/sm_89: an experiment aliasing this to std::complex
// made the whole regenerated ZA3_gpu kernel return EXACTLY 0.0 on device while the identical
// __host__ call computed correct values. Do not retry; the traces header's complex type is
// instead overridden via NT_TRACE_COMPLEX below.
#if defined(__CUDACC__) || defined(__CUDA__)
template <class T> using complex = ::cuda::std::complex<T>;
// The traces header types its complex returns via NT_TRACE_COMPLEX (default std::complex);
// point it at the SAME device-safe type so the COEN assembly's complex coefficients and the
// trace values share one complex type (mixing the two has no operator+ and, worse, std::complex
// arithmetic is the device-zeroing path this alias exists to avoid).
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX ::cuda::std::complex<double>
#endif
#else
template <class T> using complex = std::complex<T>;
#endif

// ---- math names used by generated kernels ----------------------------------
namespace compute {
template <int n> NT_DH inline double powr(double x) {
  if constexpr (n < 0) {
    double r = 1.0;
    for (int i = 0; i < -n; ++i) r *= x;
    return 1.0 / r;
  } else {
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
  }
}
// std::sqrt/pow/fma on double are the CUDA-decorated ::sqrt/::pow/::fma re-exported by
// <cmath>, so they are device-callable as-is under nvcc/clang-cuda.
using std::fma;
using std::pow;
using std::sqrt;
// complex fused multiply-add (no std::fma for complex); on the DiFfRG::complex alias so it
// binds cuda::std::complex in CUDA TUs and std::complex in host TUs.
NT_DH inline DiFfRG::complex<double> fma(const DiFfRG::complex<double> &a,
                                         const DiFfRG::complex<double> &b,
                                         const DiFfRG::complex<double> &c) {
  return a * b + c;
}
} // namespace compute

// ---- type stubs the generated kernel signatures mention --------------------
struct GPU_memory {};
template <class T> struct LogarithmicCoordinates1D {};

} // namespace DiFfRG

// ---- fixed dressing functions, shared by both paths ------------------------
// Formulas verbatim from tests/refshim/shim.hpp.
namespace shimdress {
NT_DH inline double ZA(double x) { return 1.0 / (1.0 + 0.5 * x * x) + 0.7; }
NT_DH inline double Zc(double x) { return 1.0 / (1.0 + 0.3 * x) + 0.9; }
NT_DH inline double ZAcbc(double x) { return 1.2 + 0.4 * std::sin(0.5 * x); }
NT_DH inline double ZA3(double x) { return 1.1 + 0.2 * std::cos(0.3 * x); }
NT_DH inline double ZA4(double x) { return 1.05 + 0.1 * std::sin(0.7 * x); }
NT_DH inline double dtZA(double x) { return -0.2 / (1.0 + x * x); }
NT_DH inline double dtZc(double x) { return -0.15 / (1.0 + 0.2 * x); }
NT_DH inline double ZAqbq1(double x) { return 1.15 + 0.25 * std::cos(0.4 * x); }
NT_DH inline double Zq(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }
NT_DH inline double Mq(double x) { return 0.3 / (1.0 + 0.1 * x) + 0.05; }
NT_DH inline double dtZq(double x) { return -0.18 / (1.0 + 0.25 * x); }
} // namespace shimdress

namespace DiFfRG {

// A 1-D "interpolator" stub dispatching on a runtime id (device-safe; see header comment).
// Aliased to the generated SplineInterpolator1D type name below.
struct Fn {
  int id = 0;
  NT_DH double operator()(double x) const {
    switch (id) {
    case 0: return shimdress::ZA3(x);
    case 1: return shimdress::ZAcbc(x);
    case 2: return shimdress::ZA4(x);
    case 3: return shimdress::dtZc(x);
    case 4: return shimdress::Zc(x);
    case 5: return shimdress::dtZA(x);
    case 6: return shimdress::ZA(x);
    case 7: return shimdress::ZAqbq1(x);
    case 8: return shimdress::Zq(x);
    case 9: return shimdress::Mq(x);
    case 10: return shimdress::dtZq(x);
    }
    return 0.0;
  }
};
template <class A, class B, class C> using SplineInterpolator1D = Fn;

// ---- a fixed analytic regulator (same on both paths) -----------------------
// The REG type parameter of the still-templated oracles; forwards to the free functions
// above. The `::` is load-bearing — an unqualified RB here would recurse into the member.
struct ShimRegulator {
  NT_DH static double RB(double k2, double p2) { return ::RB(k2, p2); }
  NT_DH static double RBdot(double k2, double p2) { return ::RBdot(k2, p2); }
  NT_DH static double RF(double k2, double p2) { return ::RF(k2, p2); }
  NT_DH static double RFdot(double k2, double p2) { return ::RFdot(k2, p2); }
  NT_DH static double dq2RB(double k2, double p2) { return ::dq2RB(k2, p2); }
  NT_DH static double dq2RF(double k2, double p2) { return ::dq2RF(k2, p2); }
};

} // namespace DiFfRG

// Bundle of interpolators, member names matching the kernel parameter names. Ids must match
// the Fn::operator() switch above.
struct DressingSet {
  DiFfRG::Fn ZA3{0};
  DiFfRG::Fn ZAcbc{1};
  DiFfRG::Fn ZA4{2};
  DiFfRG::Fn dtZc{3};
  DiFfRG::Fn Zc{4};
  DiFfRG::Fn dtZA{5};
  DiFfRG::Fn ZA{6};
  DiFfRG::Fn ZAqbq1{7};
  DiFfRG::Fn Zq{8};
  DiFfRG::Fn Mq{9};
  DiFfRG::Fn dtZq{10};
};
