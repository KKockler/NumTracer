// G1 precision-ladder probe types. NT_REAL_MODE selects the scalar the retyped kernel runs in:
//   1 = wrapped fp32 (float; double literals demote on construction — no 'f' suffixing needed)
//   2 = double-single (two-float compensated, ~2^-44 effective mantissa)
// The wrapper approach exists because the emitted kernels are full of double literals; bare
// `float` arithmetic would silently re-promote every mixed expression back to double.
#pragma once
#include <cmath>

#ifndef NT_REAL_MODE
#define NT_REAL_MODE 1
#endif

#define NTR_HD __host__ __device__ inline

#if NT_REAL_MODE == 1

struct nt_real {
  float v;
  nt_real() = default;
  NTR_HD nt_real(double d) : v(static_cast<float>(d)) {}
  NTR_HD explicit operator double() const { return static_cast<double>(v); }
  NTR_HD double real() const { return static_cast<double>(v); }
  NTR_HD double imag() const { return 0.0; }
};
NTR_HD nt_real operator+(nt_real a, nt_real b) { return nt_real(double(a.v + b.v)); }
NTR_HD nt_real operator-(nt_real a, nt_real b) { return nt_real(double(a.v - b.v)); }
NTR_HD nt_real operator*(nt_real a, nt_real b) { return nt_real(double(a.v * b.v)); }
NTR_HD nt_real operator/(nt_real a, nt_real b) { return nt_real(double(a.v / b.v)); }
NTR_HD nt_real operator-(nt_real a) { return nt_real(double(-a.v)); }
NTR_HD nt_real &operator+=(nt_real &a, nt_real b) { a.v += b.v; return a; }
NTR_HD nt_real &operator*=(nt_real &a, nt_real b) { a.v *= b.v; return a; }
NTR_HD nt_real sqrt(nt_real a) { return nt_real(double(sqrtf(a.v))); }
NTR_HD nt_real cos(nt_real a) { return nt_real(double(cosf(a.v))); }
NTR_HD nt_real sin(nt_real a) { return nt_real(double(sinf(a.v))); }
NTR_HD nt_real pow(nt_real a, double e) { return nt_real(double(powf(a.v, float(e)))); }
NTR_HD nt_real fma(nt_real a, nt_real b, nt_real c) { return nt_real(double(fmaf(a.v, b.v, c.v))); }
NTR_HD double nt_to_double(nt_real a) { return double(a.v); }

#elif NT_REAL_MODE == 2

// Double-single ("float-float"): Dekker/Knuth error-free transforms on the fp32 pipe.
struct nt_real {
  float hi, lo;
  nt_real() = default;
  NTR_HD nt_real(double d) : hi(static_cast<float>(d)), lo(static_cast<float>(d - double(static_cast<float>(d)))) {}
  NTR_HD nt_real(float h, float l) : hi(h), lo(l) {}
  NTR_HD explicit operator double() const { return double(hi) + double(lo); }
  NTR_HD double real() const { return double(hi) + double(lo); }
  NTR_HD double imag() const { return 0.0; }
};
NTR_HD nt_real nt_two_sum(float a, float b)
{
  const float s = a + b, bb = s - a;
  return nt_real(s, (a - (s - bb)) + (b - bb));
}
NTR_HD nt_real nt_two_prod(float a, float b)
{
  const float p = a * b;
  return nt_real(p, fmaf(a, b, -p));
}
NTR_HD nt_real operator+(nt_real a, nt_real b)
{
  nt_real s = nt_two_sum(a.hi, b.hi);
  s.lo += a.lo + b.lo;
  const float h = s.hi + s.lo;
  return nt_real(h, s.lo - (h - s.hi));
}
NTR_HD nt_real operator-(nt_real a) { return nt_real(-a.hi, -a.lo); }
NTR_HD nt_real operator-(nt_real a, nt_real b) { return a + (-b); }
NTR_HD nt_real operator*(nt_real a, nt_real b)
{
  nt_real p = nt_two_prod(a.hi, b.hi);
  p.lo += a.hi * b.lo + a.lo * b.hi;
  const float h = p.hi + p.lo;
  return nt_real(h, p.lo - (h - p.hi));
}
NTR_HD nt_real operator/(nt_real a, nt_real b)
{
  const float q1 = a.hi / b.hi;
  nt_real r = a - nt_real(q1, 0.f) * b;
  const float q2 = (r.hi + r.lo) / b.hi;
  const float h = q1 + q2;
  return nt_real(h, q2 - (h - q1));
}
NTR_HD nt_real &operator+=(nt_real &a, nt_real b) { a = a + b; return a; }
NTR_HD nt_real &operator*=(nt_real &a, nt_real b) { a = a * b; return a; }
NTR_HD nt_real sqrt(nt_real a)
{
  const float s = sqrtf(a.hi);
  if (s == 0.f) return nt_real(0.f, 0.f);
  nt_real t = nt_two_prod(s, s);
  const float e = (a.hi - t.hi - t.lo + a.lo) * 0.5f / s;
  const float h = s + e;
  return nt_real(h, e - (h - s));
}
// transcendentals: fp32 argument, correction ignored (they enter via smooth dressings only)
NTR_HD nt_real cos(nt_real a) { return nt_real(double(cosf(a.hi))); }
NTR_HD nt_real sin(nt_real a) { return nt_real(double(sinf(a.hi))); }
NTR_HD nt_real pow(nt_real a, double e) { return nt_real(double(powf(a.hi, float(e)))); }
NTR_HD nt_real fma(nt_real a, nt_real b, nt_real c) { return a * b + c; }
NTR_HD double nt_to_double(nt_real a) { return double(a.hi) + double(a.lo); }

#endif

// minimal complex carrier for the imag-only traces (ntRe/ntIm shims below). The consumer's
// own ntRe/ntIm templates call .real()/.imag(), so both types provide them as members too.
struct nt_complex {
  nt_real re, im;
  NTR_HD double real() const { return nt_to_double(re); }
  NTR_HD double imag() const { return nt_to_double(im); }
};
NTR_HD double nt_real_member_real(nt_real x) { return nt_to_double(x); }
NTR_HD double ntRe(nt_complex z) { return nt_to_double(z.re); }
NTR_HD double ntIm(nt_complex z) { return nt_to_double(z.im); }
NTR_HD double ntRe(nt_real x) { return nt_to_double(x); }
NTR_HD double ntIm(nt_real) { return 0.0; }
